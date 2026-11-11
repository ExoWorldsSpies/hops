"""Re-create tests/golden/reduction.json from the current code.

Only for intended changes of the reduction results; the reference must otherwise stay as it is.

    python tests/make_golden.py
"""

import os
import sys
import json
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# conftest sets a temporary home and puts the repository on the path, as for the tests
from conftest import TESTS  # noqa: E402
import reduction_run  # noqa: E402


def main():
    result = reduction_run.run_reduction(tempfile.mkdtemp(prefix='hops_golden_'))
    os.makedirs(os.path.join(TESTS, 'golden'), exist_ok=True)
    path = os.path.join(TESTS, 'golden', 'reduction.json')
    with open(path, 'w') as f:
        json.dump(result, f, separators=(',', ':'))
    print('{0} frames written to {1}'.format(len(result), path))


if __name__ == '__main__':
    main()
