# EECE 306 Toolbox, team starter

Rename this folder to `EECE306-Toolbox` and keep it under version control or in
a dated backup. Everything you build this semester lives here.

## Layout

```
EECE306-Toolbox/
├── +em/          the toolbox package, one subfolder per module
├── tests/        your test files, test\\\\\\\_lab01.m onward
├── runTests.m    the test runner, provided, do not modify
└── README.md     this file, replace with your own documentation
```

## First steps

1. Put the folder containing `+em` on the Octave path with `addpath`. Do **not** add `+em` itself.
2. Delete each `CONTENTS.m` placeholder as you add real functions to that module.
3. Run `runTests` from the toolbox root. With no functions written yet it will
report a failure, which is the correct starting state.

## What to put in this README by the end of the semester

Team name and members, a one line description of every module, and a short
example showing somebody else how to compute a field with your library.
The engineering quality portion of the grade is assessed partly by whether
another team could use your toolbox from this file alone.

