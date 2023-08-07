"""

.. module:: conftest.py
   :synopsis: pytest configuration
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import pytest
import os
# import sys
collect_ignore = ["__defunct__"]
# if sys.version_info[0] > 2:
#     collect_ignore.append("pkg/module_py2.py")

@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    """Causes each test to run in the directory in which the module is found **or** a subdirectory with the same base name as the module

    For example, say you have <packagename>/tests/unit/test_foo.py.  If the directory
    <packagename>/tests/unit/test_foo/ exists, the test will run in that directory; if not, it will run in <packagename>/tests/unit/

    :param request: pytest request
    :type request: pytest request
    :param monkeypatch: pytest monkeypatch
    :type monkeypatch: pytest.Monkeypatch
    """
    module_bn=request.fspath.basename
    module_name,_=os.path.splitext(module_bn)
    # print('module name',module_name)
    test_dir=request.fspath.dirname
    subdir=os.path.join(test_dir,module_name)
    if os.path.isdir(subdir):
        monkeypatch.chdir(subdir)
    else:
        monkeypatch.chdir(test_dir)
