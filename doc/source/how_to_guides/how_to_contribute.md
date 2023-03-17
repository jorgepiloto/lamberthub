# How to contribute

Main way of contributing to the library is through its usage and maintenance.
When doing so, it is possible for new issues or even bugs to appear. Those are
listed in the [official issues
board](https://github.com/jorgepiloto/lamberthub/issues).


## Developer installation

Previous situation requires for user to have `lamberthub` installed in developer
mode, so modifications can be applied to original source code. To install as
*developer*, please follow these steps:

Start by making a fork of the [official
repository](https://github.com/jorgepiloto/lamberthub) into your local machine
by running:

```bash
git clone https://github.com/your_github_username/lamberthub && cd lamberthub
```

Next, you must create a Python virtual environment, so the library dependencies
do not interfere with your system environment:

```bash
python -m venv .venv && source .venv/bin/activate
```

Check your current Python binary is the one hosted in the virtual environment by
running:

```bash
which python
```

which must retrieve the absolute path of the previously created virtual
environment:

```bash
/lamberthub/.venv/bin/python
```

Now, update `pip`:

```bash
python -m pip install -U pip
```

Finally, install `lamberthub` in development mode:

```bash
python -m pip install --editable .
```

which tells your Python environment that the `lamberthub` package is located in
this actual folder. **Now, any change you apply to source code affects the
behavior of the library**.
