"""

.. module:: test_configuration
   :synopsis: unit tests for htpolynet.core.configuration

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import os
import json
import tempfile
import shutil
import logging
logger = logging.getLogger(__name__)

import yaml
from htpolynet.core.configuration import Configuration


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_yaml(path, data):
    with open(path, 'w') as f:
        yaml.dump(data, f)


def _write_json(path, data):
    with open(path, 'w') as f:
        json.dump(data, f)


MINIMAL_DATA = {
    'Title': 'Test system',
    'constituents': {
        'MOL': {'count': 10},
    },
}

# Every key here is real.  This fixture used to carry invented ones --
# 'mdrun_opts', 'CURE.max_conversion', 'GAFF.minimize_molecules',
# 'precure.equilibration' -- which the old parser discarded in silence, so the
# tests passed while asserting nothing about a real configuration.  Validation
# rejects them now, which is how they were found.
FULL_DATA = {
    'Title': 'Full test system',
    'ncpu': 4,
    'constituents': {
        'EPA': {'count': 50, 'stereocenters': []},
        'EPI': {'count': 25},
    },
    'reactions': [
        {
            'name': 'rxn1',
            'stage': 'cure',
            'reactants': {1: 'EPA', 2: 'EPI'},
            'product': 'EPA~C1-N1~EPI',
            'atoms': {
                'A': {'reactant': 1, 'resid': 1, 'atom': 'C1', 'z': 1},
                'B': {'reactant': 2, 'resid': 1, 'atom': 'N1', 'z': 1},
            },
            'bonds': [{'atoms': ['A', 'B'], 'order': 1}],
        },
    ],
    'gromacs': {'gmx': 'gmx_mpi', 'mdrun_options': {'ntmpi': 1}},
    'densification': {'aspect_ratio': [1.0, 1.0, 1.0], 'scale': 0.4},
    'precure': {'preequilibration': {'ensemble': 'npt', 'ps': 50}},
    'postcure': {'postequilibration': {'ensemble': 'nvt', 'ps': 50}},
    'CURE': {'controls': {'desired_conversion': 0.85, 'max_iterations': 10}},
    'GAFF': {'resolve_type_discrepancies': [
        {'typename': 'dihedraltypes', 'funcidx': 4, 'rule': 'stiffest'}]},
    'ambertools': {'charge_method': 'bcc'},
    'resolve_type_discrepancies': [
        {'typename': 'dihedraltypes', 'funcidx': 4, 'rule': 'stiffest'}],
}


class TestConfigurationDefaults(unittest.TestCase):
    """__init__ sets sensible defaults before any file is read."""

    def setUp(self):
        self.c = Configuration()

    def test_cfgfile_empty(self):
        self.assertEqual(self.c.cfgfile, '')

    def test_title_empty(self):
        self.assertEqual(self.c.title, '')

    def test_ncpu_defaults_to_cpu_count(self):
        self.assertEqual(self.c.ncpu, os.cpu_count())

    def test_constituents_empty_dict(self):
        self.assertIsInstance(self.c.constituents, dict)
        self.assertEqual(len(self.c.constituents), 0)

    def test_reaction_specs_empty_list(self):
        self.assertIsInstance(self.c.reaction_specs, list)
        self.assertEqual(len(self.c.reaction_specs), 0)

    def test_initial_composition_empty_list(self):
        self.assertIsInstance(self.c.initial_composition, list)
        self.assertEqual(len(self.c.initial_composition), 0)

    def test_gromacs_empty_dict(self):
        self.assertIsInstance(self.c.gromacs, dict)

    def test_all_phase_dicts_empty(self):
        for attr in ('densification', 'precure', 'postcure', 'cure', 'gaff', 'ambertools'):
            self.assertIsInstance(getattr(self.c, attr), dict, f'{attr} should be a dict')
            self.assertEqual(len(getattr(self.c, attr)), 0, f'{attr} should be empty')

    def test_resolve_type_discrepancies_empty_list(self):
        self.assertIsInstance(self.c.resolve_type_discrepancies, list)

    def test_basedict_empty(self):
        self.assertIsInstance(self.c.basedict, dict)
        self.assertEqual(len(self.c.basedict), 0)


class TestConfigurationReadYAML(unittest.TestCase):
    """Configuration.read() with YAML files."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _yaml_path(self, name='cfg.yaml'):
        return os.path.join(self.tmpdir, name)

    # --- basic parsing ---

    def test_read_yaml_returns_configuration(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertIsInstance(c, Configuration)

    def test_cfgfile_set(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.cfgfile, p)

    def test_title_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.title, 'Test system')

    def test_yml_extension_accepted(self):
        p = self._yaml_path('cfg.yml')
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertIsInstance(c, Configuration)

    # --- constituents ---

    def test_constituents_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertIn('MOL', c.constituents)
        self.assertEqual(c.constituents['MOL']['count'], 10)

    def test_multiple_constituents(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(set(c.constituents.keys()), {'EPA', 'EPI'})

    # --- initial_composition derived from constituents ---

    def test_initial_composition_length(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(len(c.initial_composition), 2)

    def test_initial_composition_keys(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        entry = c.initial_composition[0]
        self.assertIn('molecule', entry)
        self.assertIn('count', entry)

    def test_initial_composition_molecule_name(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.initial_composition[0]['molecule'], 'MOL')

    def test_initial_composition_count(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.initial_composition[0]['count'], 10)

    def test_initial_composition_count_missing_defaults_zero(self):
        data = {'constituents': {'X': {}}}
        p = self._yaml_path()
        _write_yaml(p, data)
        c = Configuration.read(p)
        self.assertEqual(c.initial_composition[0]['count'], 0)

    # --- reactions ---

    def test_reaction_specs_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(len(c.reaction_specs), 1)
        self.assertEqual(c.reaction_specs[0]['name'], 'rxn1')

    def test_reaction_specs_absent_defaults_empty(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.reaction_specs, [])

    # --- optional sections ---

    def test_gromacs_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.gromacs.get('gmx'), 'gmx_mpi')

    def test_gromacs_absent_gets_schema_defaults(self):
        # the old parser left an omitted section empty and let Runtime fill it
        # in later; the schema is now the single place those defaults live
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.gromacs['gmx'], 'gmx')
        self.assertEqual(c.gromacs['gmx_options'], '-quiet -nobackup')

    def test_densification_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertAlmostEqual(c.densification.get('scale'), 0.4)

    def test_precure_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertIn('preequilibration', c.precure)
        self.assertEqual(c.precure['preequilibration']['ps'], 50)

    def test_postcure_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertIn('postequilibration', c.postcure)
        self.assertEqual(c.postcure['postequilibration']['ensemble'], 'nvt')

    def test_cure_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertAlmostEqual(c.cure['controls']['desired_conversion'], 0.85)
        self.assertEqual(c.cure['controls']['max_iterations'], 10)

    def test_cure_unset_controls_get_their_documented_defaults(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertAlmostEqual(c.cure['controls']['search_radius'], 0.5)
        self.assertAlmostEqual(c.cure['relax']['increment'], 0.08)

    def test_gaff_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.gaff['resolve_type_discrepancies'][0]['rule'], 'stiffest')

    def test_ambertools_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.ambertools.get('charge_method'), 'bcc')

    def test_resolve_type_discrepancies_parsed(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(len(c.resolve_type_discrepancies), 1)
        self.assertEqual(c.resolve_type_discrepancies[0]['typename'], 'dihedraltypes')

    # --- ncpu ---

    def test_ncpu_from_file(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.ncpu, 4)

    def test_ncpu_absent_defaults_to_cpu_count(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.ncpu, os.cpu_count())

    # --- basedict preserved ---

    def test_basedict_preserved(self):
        p = self._yaml_path()
        _write_yaml(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.basedict['Title'], 'Full test system')

    def test_basedict_is_raw_dict(self):
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertIsInstance(c.basedict, dict)

    # --- missing optional sections default to empty ---

    def test_missing_sections_get_documented_defaults(self):
        """An omitted section is filled from the schema, not left empty."""
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        for attr in ('gromacs', 'densification', 'precure', 'postcure', 'cure', 'ambertools'):
            self.assertTrue(getattr(c, attr), f'{attr} should be populated from the schema')
        self.assertAlmostEqual(c.densification['initial_density'], 200.0)
        self.assertEqual(c.ambertools['charge_method'], 'gas')

    def test_presence_sensitive_keys_are_not_injected(self):
        """Keys whose presence changes behavior must stay absent.

        runtime.py asserts that exactly one of initial_density and
        initial_boxsize is *present*, so filling in an empty initial_boxsize
        would fail that assertion on every build.
        """
        p = self._yaml_path()
        _write_yaml(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertNotIn('initial_boxsize', c.densification)
        self.assertNotIn('reactive_atoms', c.constituents['MOL'])
        self.assertNotIn('rename_atoms', c.constituents['MOL'])
        self.assertEqual(c.resolve_type_discrepancies, [])


class TestConfigurationReadJSON(unittest.TestCase):
    """Configuration.read() with JSON files."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _json_path(self):
        return os.path.join(self.tmpdir, 'cfg.json')

    def test_read_json_returns_configuration(self):
        p = self._json_path()
        _write_json(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertIsInstance(c, Configuration)

    def test_json_title_parsed(self):
        p = self._json_path()
        _write_json(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.title, 'Test system')

    def test_json_constituents_parsed(self):
        p = self._json_path()
        _write_json(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertIn('EPA', c.constituents)

    def test_json_initial_composition_correct(self):
        p = self._json_path()
        _write_json(p, MINIMAL_DATA)
        c = Configuration.read(p)
        self.assertEqual(len(c.initial_composition), 1)
        self.assertEqual(c.initial_composition[0]['molecule'], 'MOL')

    def test_json_reactions_parsed(self):
        p = self._json_path()
        _write_json(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(len(c.reaction_specs), 1)

    def test_json_gromacs_parsed(self):
        p = self._json_path()
        _write_json(p, FULL_DATA)
        c = Configuration.read(p)
        self.assertEqual(c.gromacs.get('gmx'), 'gmx_mpi')


class TestConfigurationErrors(unittest.TestCase):
    """Error handling in Configuration.read()."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_unknown_extension_raises(self):
        p = os.path.join(self.tmpdir, 'cfg.txt')
        _write_yaml(p, MINIMAL_DATA)
        with self.assertRaises(Exception):
            Configuration.read(p)

    def test_missing_file_raises(self):
        p = os.path.join(self.tmpdir, 'nonexistent.yaml')
        with self.assertRaises(FileNotFoundError):
            Configuration.read(p)

    def test_bad_yaml_raises(self):
        p = os.path.join(self.tmpdir, 'bad.yaml')
        with open(p, 'w') as f:
            f.write(': : : invalid yaml :::\n')
        with self.assertRaises(Exception):
            Configuration.read(p)

    def test_bad_json_raises(self):
        p = os.path.join(self.tmpdir, 'bad.json')
        with open(p, 'w') as f:
            f.write('{not valid json')
        with self.assertRaises(Exception):
            Configuration.read(p)


class TestConfigurationTitleFallback(unittest.TestCase):
    """Title key is optional; default text is provided when absent."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def test_title_absent_gives_default(self):
        data = {'constituents': {'A': {'count': 1}}}
        p = os.path.join(self.tmpdir, 'cfg.yaml')
        _write_yaml(p, data)
        c = Configuration.read(p)
        self.assertNotEqual(c.title, '')
        self.assertIsInstance(c.title, str)

    def test_title_present_overrides_default(self):
        data = {'Title': 'My polymer', 'constituents': {}}
        p = os.path.join(self.tmpdir, 'cfg.yaml')
        _write_yaml(p, data)
        c = Configuration.read(p)
        self.assertEqual(c.title, 'My polymer')


class TestConfigurationInitialComposition(unittest.TestCase):
    """initial_composition is consistently derived from constituents."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir)

    def _read(self, data):
        p = os.path.join(self.tmpdir, 'cfg.yaml')
        _write_yaml(p, data)
        return Configuration.read(p)

    def test_empty_constituents_gives_empty_composition(self):
        c = self._read({'constituents': {}})
        self.assertEqual(c.initial_composition, [])

    def test_composition_entries_match_constituent_names(self):
        data = {'constituents': {'A': {'count': 5}, 'B': {'count': 3}}}
        c = self._read(data)
        names = {e['molecule'] for e in c.initial_composition}
        self.assertEqual(names, {'A', 'B'})

    def test_composition_counts_match_constituents(self):
        data = {'constituents': {'X': {'count': 7}}}
        c = self._read(data)
        self.assertEqual(c.initial_composition[0]['count'], 7)

    def test_composition_and_constituents_stay_in_sync(self):
        """Every constituent name must appear in initial_composition."""
        data = {
            'constituents': {
                'MON': {'count': 100},
                'CRS': {'count': 10},
                'CAP': {'count': 20},
            }
        }
        c = self._read(data)
        comp_names = {e['molecule'] for e in c.initial_composition}
        self.assertEqual(comp_names, set(c.constituents.keys()))


class TestConfigurationValidation(unittest.TestCase):
    """What the schema rejects that the old parser accepted in silence."""

    def setUp(self):
        self.d = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _read(self, mutate):
        import copy
        data = copy.deepcopy(FULL_DATA)
        mutate(data)
        p = os.path.join(self.d, 'cfg.yaml')
        _write_yaml(p, data)
        return Configuration.read(p)

    def test_misspelled_control_is_rejected(self):
        # the costly case: desired_conversion silently ignored is a half-cure
        with self.assertRaises(Exception) as cm:
            self._read(lambda d: d['CURE']['controls'].__setitem__(
                'desired_converson', d['CURE']['controls'].pop('desired_conversion')))
        self.assertIn('desired_converson', str(cm.exception))

    def test_unknown_top_level_section_is_rejected(self):
        with self.assertRaises(Exception) as cm:
            self._read(lambda d: d.__setitem__('postcure_repairs', []))
        self.assertIn('postcure_repairs', str(cm.exception))

    def test_bad_ensemble_is_rejected(self):
        with self.assertRaises(Exception) as cm:
            self._read(lambda d: d['precure']['preequilibration'].__setitem__('ensemble', 'nvE'))
        self.assertIn('nvE', str(cm.exception))

    def test_wrong_type_is_rejected(self):
        with self.assertRaises(Exception) as cm:
            self._read(lambda d: d['constituents']['EPA'].__setitem__('count', 'fifty'))
        self.assertIn('count', str(cm.exception))

    def test_the_error_names_the_config_file(self):
        with self.assertRaises(Exception) as cm:
            self._read(lambda d: d['gromacs'].__setitem__('mdrun_opts', '-ntmpi 1'))
        self.assertIn('cfg.yaml', str(cm.exception))


class TestSchemaMatchesPythonDefaults(unittest.TestCase):
    """The schema and the Python default dicts must not drift apart.

    They are two statements of the same defaults, and this suite has already
    caught six cases where a documented default disagreed with the code.  While
    both exist, something has to compare them; the moment relax.increment was
    fixed in Python only, this test is what noticed.
    """

    def _schema(self):
        from htpolynet.core.configuration import schema_path
        with schema_path() as p:
            with open(p) as f:
                return yaml.safe_load(f)

    def _leaf_defaults(self, attrs, prefix=''):
        out = {}
        for a in attrs or []:
            name = f"{prefix}{a['name']}"
            kids = a.get('attributes')
            if kids:
                out.update(self._leaf_defaults(kids, f'{name}.'))
            elif 'default' in a:
                out[name] = a['default']
        return out

    def _section(self, section):
        top = {a['name']: a for a in self._schema()['attributes']}
        return self._leaf_defaults(top[section].get('attributes'))

    def test_cure_defaults_agree(self):
        from htpolynet.cure.curecontroller import CureController
        schema = self._section('CURE')
        for block, d in CureController.curedict_defaults.items():
            for k, v in d.items():
                key = f'{block}.{k}'
                if key not in schema or isinstance(v, list):
                    continue          # lists compared separately; ncpu is host-dependent
                if k == 'ncpu':
                    continue
                self.assertEqual(schema[key], v, f'{key}: schema {schema[key]!r} vs code {v!r}')

    def test_gromacs_and_ambertools_defaults_agree(self):
        from htpolynet.core.runtime import Runtime
        for section, cfgkey in (('gromacs', 'gromacs'), ('ambertools', 'ambertools')):
            schema = self._section(section)
            for k, v in Runtime.runtime_defaults[cfgkey].items():
                if k in schema:
                    self.assertEqual(schema[k], v, f'{section}.{k}')

    def test_densification_initial_density_agrees(self):
        from htpolynet.core.runtime import Runtime
        self.assertEqual(self._section('densification')['initial_density'],
                         Runtime.runtime_defaults['densification']['initial_density'])
