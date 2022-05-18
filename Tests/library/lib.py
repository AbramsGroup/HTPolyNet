import importlib.resources
import os
package='Library'
x=importlib.resources.contents(package)
for g in x:
    print(g)
with importlib.resources.path(package,'__init__.py') as f:
    print(os.path.abspath(f))

