import os
import sys
import tempfile

TESTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(TESTS))
sys.path.insert(0, TESTS)

# HOPSLog reads and writes ~/hops_log_*.yaml: use a temporary home, so the tests start from the
# package defaults and do not touch the user's settings
os.environ['HOME'] = tempfile.mkdtemp(prefix='hops_test_home_')
