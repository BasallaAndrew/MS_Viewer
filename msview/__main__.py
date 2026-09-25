"""Entry point: msview [file]. With no file, a file picker opens at startup."""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="msview", description="Browse Thermo .raw and mzML LC-MS files.")
    parser.add_argument("file", nargs="?", help="a .raw or .mzML file to open")
    parser.add_argument("--no-picker", action="store_true", help="don't show the file picker at startup")
    args = parser.parse_args(argv)

    from PySide6 import QtCore, QtWidgets

    from msview.app import MainWindow

    app = QtWidgets.QApplication(sys.argv[:1])
    app.setOrganizationName("msview")
    app.setApplicationName("msview")
    window = MainWindow()
    window.show()
    window.raise_()
    if args.file:
        QtCore.QTimer.singleShot(0, lambda: window.open_file(args.file))
    elif not args.no_picker:
        QtCore.QTimer.singleShot(0, window.choose_file)
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
