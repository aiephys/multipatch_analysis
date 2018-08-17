"""
Runs all pipeline analysis stages on a daily schedule. 

Why not use cron like a normal human being?
Because I like having the script run in the console where I cam monitor and
debug problems more easily.
"""

import os, sys, time, argparse
from datetime import datetime, timedelta

def delay(hour=2):
    """Sleep until *hour*"""
    now = datetime.now()
    tomorrow = now + timedelta(days=1)
    next_run = datetime(tomorrow.year, tomorrow.month, tomorrow.day, 2, 0)
    delay = (next_run - now).total_seconds()

    print("Sleeping %d seconds until %s.." % (delay, next_run))
    time.sleep(delay)


def run(cmd, msg):
    print("======================================================================================")
    print("    " + msg)
    print("======================================================================================")
    os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run all analysis pipeline stages to import / analyze new data on a schedule.")
    parser.add_argument('--now', default=False, action='store_true', help="Run once immediately before starting scheduled updates.")

    args = parser.parse_args(sys.argv[1:])

    if not args.now:
        delay()

    while True:
        run('python util/sync_rigs_to_server.py', 'sync raw data')
        run('python util/import_to_database.py', 'import to DB')
        run('python util/update_morphology.py', 'update morphology')
        run('python util/analyze_pulse_response_strength.py', 'pulse response strength')
        run('python util/analyze_connection_strength.py', 'connection strength')
        run('python util/database.py --vacuum', 'vacuum')

        delay()


