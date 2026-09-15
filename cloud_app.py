"""Community Cloud entry point: bootstrap the same pipeline on a fresh deployment."""

import runpy

from analysis import run_analysis
from load_data import DB_PATH, ROOT, load_data


def initialize():
    if not DB_PATH.exists():
        load_data()
        run_analysis()


initialize()
runpy.run_path(str(ROOT / "app.py"), run_name="__main__")
