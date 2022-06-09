.. _tutorial_pms_monomer:

monomer
=======

In this section, we describe how the input :download:`MST.mol2 <MST.mol2>` that specifies the 4-methylstyrene monomer is generated.  Since this represents an instance where a new system is being generated, let's begin by creating an empty directory and then populating with a "molecule library":abbr:

.. code-block:: console

    $ cd 
    $ mkdir my_pms_build
    $ cd my_pms_build
    $ mkdir lib
    $ mkdir lib/inputs
    $ mkdir lib/parameterized
    $ cd lib/inputs

Now we can generate the required ``*.mol2`` file.