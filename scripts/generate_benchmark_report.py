"""Generate the README performance comparison from benchmark JSON reports."""

from __future__ import annotations

import argparse
from datetime import UTC, datetime
import json
from pathlib import Path
import re

START_MARKER = "<!-- performance-comparison:start -->"
END_MARKER = "<!-- performance-comparison:end -->"
SECTION_TITLE = "## Performance comparison"

GROUP_TITLES = {
    "zero-rev-nominal": "Zero-revolution nominal cases",
    "zero-rev-general": "Zero-revolution general cases",
    "zero-rev-hyperbolic": "Zero-revolution hyperbolic cases",
    "multi-rev-m1": "One-revolution branch cases",
    "multi-rev-hard": "Near-minimum-energy multi-revolution cases",
    "robust-zero-rev": "Near-tangent robustness cases",
}

CASE_TITLE_OVERRIDES = {
    "battin-book": "Battin book",
    "curtiss-book": "Curtiss book",
    "vallado-book": "Vallado book",
    "gmat-prograde": "GMAT hyperbolic",
    "gmat-retrograde": "GMAT hyperbolic",
    "der-article-i-prograde-low": "Der article I",
    "der-article-i-retrograde-high": "Der article I",
    "der-article-ii-prograde-high": "Der article II",
    "der-article-ii-retrograde-high": "Der article II",
    "prograde-high": "Der article II",
    "prograde-low": "Der article II",
    "retrograde-high": "Der article II",
    "retrograde-low": "Der article II",
    "m1-prograde-high": "Near-minimum-energy",
    "m1-prograde-low": "Near-minimum-energy",
    "m2-prograde-high": "Near-minimum-energy",
    "m2-prograde-low": "Near-minimum-energy",
    "near-tangent-prograde-high": "Near tangent",
}


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--benchmarks-dir",
        type=Path,
        default=Path(".benchmarks"),
        help="Directory containing pytest-benchmark JSON reports.",
    )
    parser.add_argument(
        "--readme",
        type=Path,
        default=Path("README.md"),
        help="README file to update.",
    )
    return parser.parse_args()


def load_reports(benchmarks_dir):
    """Load all pytest-benchmark JSON reports from a directory."""
    reports = []
    for path in sorted(benchmarks_dir.glob("*.json")):
        with path.open(encoding="utf-8") as file_handler:
            reports.append(json.load(file_handler))
    if not reports:
        raise FileNotFoundError(f"No benchmark JSON reports found in {benchmarks_dir}")
    return reports


def solver_name(benchmark):
    """Extract a solver name from a pytest-benchmark entry."""
    solver = benchmark["params"]["solver"]
    match = re.search(r"<function ([a-zA-Z0-9_]+) ", solver)
    if match:
        return match.group(1)
    return solver


def case_key(benchmark):
    """Return a stable key for one benchmark case."""
    case = benchmark["params"]["case"]
    return (
        case["name"],
        case["M"],
        case["is_prograde"],
        case["is_low_path"],
    )


def case_title(case):
    """Return a readable case title."""
    name = case["name"]
    if name in CASE_TITLE_OVERRIDES:
        return CASE_TITLE_OVERRIDES[name]
    return name.replace("-", " ").title()


def case_prograde(case):
    """Return a Yes/No label for the prograde flag."""
    return "Yes" if case["is_prograde"] else "No"


def case_path(case):
    """Return the selected Lambert path."""
    return "low" if case["is_low_path"] else "high"


def group_title(group):
    """Return a readable benchmark group title."""
    return GROUP_TITLES.get(group, group.replace("-", " ").title())


def format_microseconds(seconds):
    """Format seconds as microseconds with one decimal place."""
    return f"{seconds * 1e6:.1f}"


def sorted_benchmarks(reports):
    """Group benchmark entries by benchmark group, sorted by case then median."""
    groups = {}
    for report in reports:
        for benchmark in report["benchmarks"]:
            groups.setdefault(benchmark["group"], []).append(benchmark)

    for benchmarks in groups.values():
        benchmarks.sort(
            key=lambda item: (
                case_title(item["params"]["case"]),
                item["params"]["case"]["M"],
                item["params"]["case"]["is_prograde"],
                item["params"]["case"]["is_low_path"],
                item["stats"]["median"],
                solver_name(item),
            )
        )
    return dict(sorted(groups.items()))


def machine_summary(report):
    """Build a compact machine summary from a pytest-benchmark report."""
    machine_info = report.get("machine_info", {})
    cpu_info = machine_info.get("cpu", {})
    cpu = cpu_info.get("brand_raw", "unknown CPU")
    python = machine_info.get("python_version", "unknown Python")
    system = machine_info.get("system", "unknown OS")
    return f"{system}, Python {python}, {cpu}"


def _rank_within_case(benchmarks):
    """Return a dict mapping (case_key, solver_name) -> rank within that case."""
    by_case: dict[tuple, list] = {}
    for b in benchmarks:
        by_case.setdefault(case_key(b), []).append(b)

    ranks = {}
    for key, entries in by_case.items():
        sorted_entries = sorted(entries, key=lambda b: b["stats"]["median"])
        for rank, entry in enumerate(sorted_entries, start=1):
            ranks[(key, solver_name(entry))] = rank
    return ranks


def _speedup_within_case(benchmarks):
    """Return a dict mapping (case_key, solver_name) -> speedup vs slowest."""
    by_case: dict[tuple, list] = {}
    for b in benchmarks:
        by_case.setdefault(case_key(b), []).append(b)

    speedups = {}
    for key, entries in by_case.items():
        slowest = max(e["stats"]["median"] for e in entries)
        for entry in entries:
            median = entry["stats"]["median"]
            speedups[(key, solver_name(entry))] = (
                slowest / median if median > 0 else 1.0
            )
    return speedups


def _benchmarks_by_scenario(benchmarks):
    """Group benchmark entries by readable scenario title."""
    by_scenario: dict[str, list] = {}
    for benchmark in benchmarks:
        by_scenario.setdefault(case_title(benchmark["params"]["case"]), []).append(
            benchmark
        )
    return by_scenario


def build_section(reports):
    """Build the generated Markdown performance section."""
    groups = sorted_benchmarks(reports)
    commit = reports[0].get("commit_info", {}).get("id", "unknown")
    generated_at = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")

    lines = [
        SECTION_TITLE,
        "",
        START_MARKER,
        "",
        "_This section is auto-generated from CI benchmark artifacts._",
        "",
        f"- Generated: {generated_at}",
        f"- Commit: `{commit[:12]}`",
        f"- Environment: {machine_summary(reports[0])}",
        "",
        "Times are in microseconds (lower is better).",
        "**Speedup** is relative to the slowest solver for each case.",
        "All solvers are JIT-warmed before timing begins.",
        "",
    ]

    for group, benchmarks in groups.items():
        ranks = _rank_within_case(benchmarks)
        speedups = _speedup_within_case(benchmarks)

        lines.extend(
            [
                f"### {group_title(group)}",
                "",
            ]
        )
        for scenario, scenario_benchmarks in _benchmarks_by_scenario(
            benchmarks
        ).items():
            lines.extend(
                [
                    f"#### {scenario}",
                    "",
                    "| Rank | Solver | Revolutions | Prograde | Path | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |",
                    "|-----:|--------|------------:|----------|------|-----------:|----------:|---------:|--------:|-------:|",
                ]
            )
            for benchmark in scenario_benchmarks:
                case = benchmark["params"]["case"]
                name = solver_name(benchmark)
                stats = benchmark["stats"]
                key = case_key(benchmark)
                rank = ranks.get((key, name), 0)
                speedup = speedups.get((key, name), 1.0)
                lines.append(
                    f"| {rank} | "
                    f"`{name}` | "
                    f"{case['M']} | "
                    f"{case_prograde(case)} | "
                    f"{case_path(case)} | "
                    f"{format_microseconds(stats['median'])} | "
                    f"{format_microseconds(stats['mean'])} | "
                    f"{format_microseconds(stats['iqr'])} | "
                    f"{speedup:.1f}x | "
                    f"{stats['rounds']} |"
                )
            lines.append("")

    lines.extend([END_MARKER, ""])
    return "\n".join(lines)


def update_readme(readme_path, section):
    """Replace or append the generated performance section in the README."""
    content = readme_path.read_text(encoding="utf-8")
    if START_MARKER in content and END_MARKER in content:
        generated_content = section.split(f"{START_MARKER}\n\n", maxsplit=1)[1]
        generated_content = f"{START_MARKER}\n\n{generated_content}"
        pattern = re.compile(
            rf"{re.escape(START_MARKER)}.*?{re.escape(END_MARKER)}\n?",
            flags=re.DOTALL,
        )
        updated = pattern.sub(generated_content, content)
    else:
        updated = content.rstrip() + "\n\n" + section

    readme_path.write_text(updated, encoding="utf-8")


def main():
    """Generate and inject the performance comparison section."""
    args = parse_args()
    reports = load_reports(args.benchmarks_dir)
    section = build_section(reports)
    update_readme(args.readme, section)


if __name__ == "__main__":
    main()
