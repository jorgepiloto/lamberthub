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
                item["params"]["case"]["name"],
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
    """Return a dict mapping (case_name, solver_name) -> rank within that case."""
    by_case: dict[str, list] = {}
    for b in benchmarks:
        case = b["params"]["case"]["name"]
        by_case.setdefault(case, []).append(b)

    ranks = {}
    for case, entries in by_case.items():
        sorted_entries = sorted(entries, key=lambda b: b["stats"]["median"])
        for rank, entry in enumerate(sorted_entries, start=1):
            ranks[(case, solver_name(entry))] = rank
    return ranks


def _speedup_within_case(benchmarks):
    """Return a dict mapping (case_name, solver_name) -> speedup vs slowest."""
    by_case: dict[str, list] = {}
    for b in benchmarks:
        case = b["params"]["case"]["name"]
        by_case.setdefault(case, []).append(b)

    speedups = {}
    for case, entries in by_case.items():
        slowest = max(e["stats"]["median"] for e in entries)
        for entry in entries:
            median = entry["stats"]["median"]
            speedups[(case, solver_name(entry))] = (
                slowest / median if median > 0 else 1.0
            )
    return speedups


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
                f"### {group}",
                "",
                "| Rank | Case | Solver | Median (µs) | Mean (µs) | IQR (µs) | Speedup | Rounds |",
                "|-----:|------|--------|-----------:|----------:|---------:|--------:|-------:|",
            ]
        )
        for benchmark in benchmarks:
            case = benchmark["params"]["case"]["name"]
            name = solver_name(benchmark)
            stats = benchmark["stats"]
            rank = ranks.get((case, name), 0)
            speedup = speedups.get((case, name), 1.0)
            lines.append(
                f"| {rank} | "
                f"`{case}` | "
                f"`{name}` | "
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
        pattern = re.compile(
            rf"{re.escape(SECTION_TITLE)}\n\n{re.escape(START_MARKER)}.*?"
            rf"{re.escape(END_MARKER)}\n?",
            flags=re.DOTALL,
        )
        updated = pattern.sub(section, content)
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
