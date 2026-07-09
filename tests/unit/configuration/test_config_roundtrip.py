"""
Contract for configuration (de)serialization and hashing.

This suite pins the round-trip guarantees a shared MolForge config must satisfy:

- One canonical, flat/step-keyed on-disk format for JSON and YAML, covering core
  actors and plugins uniformly.
- Every entry point (``MolForge(config_path=...)``, ``PipeParams.from_json`` /
  ``from_yaml`` / ``from_dict``, and ``build_params``) reads that format.
- A config object -> file -> config object round-trip reproduces identical
  parameters and an identical hash.
- Run metadata (run_id, log_time, input_id, ...) may be recorded in a written
  config, but is ignored when a params object is built from it and excluded from
  the hash.
- Private attributes (keys starting with ``_``) never reach the file.
- The MolForge version participates in the hash and is defined in one place, so a
  hash match implies a version match.
"""

import json
from dataclasses import dataclass, field

import pytest
import yaml

import molforge
from molforge import ForgeParams, MolForge
from molforge.actors.params.base import BaseParams
from molforge.configuration.factory import PipeParams
from molforge.configuration.serialization import ConfigCodec


STEPS = ['source', 'chembl', 'curate', 'tokens', 'distributions', 'scaffold', 'split']


# Synthetic nested params for the recursive-reconstruction contract. A plugin
# author may nest a params dataclass inside another; the inner object carries a
# private attribute derived in ``__post_init__`` that stays out of the config and
# is re-derived on load rather than injected.

@dataclass
class _InnerParams(BaseParams):
    scale: float = 1.0

    def _validate_params(self) -> None:
        pass

    def _post_init_hook(self) -> None:
        self._scaled = self.scale * 10  # derived, private


@dataclass
class _OuterParams(BaseParams):
    label: str = 'base'
    inner: _InnerParams = field(default_factory=_InnerParams)

    def _validate_params(self) -> None:
        pass


@dataclass
class _OuterUnionParams(BaseParams):
    # PEP 604 ``X | None`` annotation on a nested-params field.
    inner: _InnerParams | None = None

    def _validate_params(self) -> None:
        pass


def make_params():
    """A representative config exercising core actors and a plugin."""
    return ForgeParams(
        steps=STEPS,
        source_params={'backend': 'sql', 'version': 36},
        chembl_params={'standard_type': 'IC50'},
        plugin_params={'split': {'test_ratio': 0.15, 'val_ratio': 0.1}},
    )


# --------------------------------------------------------------------------- #
# Canonical format
# --------------------------------------------------------------------------- #

class TestCanonicalFormat:

    @pytest.mark.parametrize('writer,loader', [('to_json', json.loads), ('to_yaml', yaml.safe_load)])
    def test_written_config_is_flat_step_keyed(self, tmp_path, writer, loader):
        p = make_params()
        path = tmp_path / 'config'
        getattr(p, writer)(str(path))
        data = loader(path.read_text())

        # Per-step blocks are keyed by step name for core AND plugin steps;
        # the internal *_params / plugin_params representation does not leak out.
        for step in STEPS:
            assert step in data, f"missing step block '{step}'"
        assert 'source_params' not in data
        assert 'plugin_params' not in data
        assert data['steps'] == STEPS

    def test_no_private_attrs_in_written_config(self, tmp_path):
        p = make_params()
        path = tmp_path / 'config.json'
        p.to_json(str(path))

        def assert_clean(obj):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    assert not k.startswith('_'), f"private key leaked: {k}"
                    assert_clean(v)

        assert_clean(json.loads(path.read_text()))


# --------------------------------------------------------------------------- #
# Round-trip identity
# --------------------------------------------------------------------------- #

class TestRoundTrip:

    @pytest.mark.parametrize('ext,writer', [('json', 'to_json'), ('yaml', 'to_yaml')])
    def test_object_file_object_identity(self, tmp_path, ext, writer):
        p = make_params()
        path = tmp_path / f'config.{ext}'
        getattr(p, writer)(str(path))
        p2 = MolForge(config_path=str(path)).config
        assert p2._get_config_hash() == p._get_config_hash()

    def test_plugin_params_survive_roundtrip(self, tmp_path):
        p = make_params()
        path = tmp_path / 'config.json'
        p.to_json(str(path))
        p2 = PipeParams.from_json(str(path))
        assert p2.plugin_params['split'].test_ratio == 0.15
        assert p2.plugin_params['split'].val_ratio == 0.1

    def test_flat_cli_config_loads_through_same_path(self):
        # The hand-written demo/template flat format loads via the Python API too.
        mf = MolForge(config_path='demo/jak2_ic50.yaml')
        assert 'split' in mf.config.steps and 'confs' in mf.config.steps


# --------------------------------------------------------------------------- #
# Cross-API passability (CLI <-> Python)
# --------------------------------------------------------------------------- #

class TestCrossApi:

    def test_cli_config_loads_in_python_api(self, tmp_path):
        # A config assembled by the CLI builder writes and reloads via the Python API.
        from molforge.cli import builder
        p = builder.build_params(steps=STEPS, overrides=['split.test_ratio=0.2'])
        path = tmp_path / 'cli.yaml'
        p.to_yaml(str(path))
        assert MolForge(config_path=str(path)).config._get_config_hash() == p._get_config_hash()

    def test_python_config_loads_in_cli_builder(self, tmp_path):
        # A config written by the Python API reloads via the CLI builder.
        from molforge.cli import builder
        p = make_params()
        path = tmp_path / 'py.yaml'
        p.to_yaml(str(path))
        assert builder.build_params(config_path=str(path))._get_config_hash() == p._get_config_hash()


# --------------------------------------------------------------------------- #
# Run metadata
# --------------------------------------------------------------------------- #

class TestRunMetadata:

    def test_run_metadata_ignored_on_load_and_hash(self, tmp_path):
        p = make_params()
        clean_path = tmp_path / 'clean.json'
        p.to_json(str(clean_path))

        # A run records provenance into its written config; loading it back must
        # ignore that metadata and reproduce the pristine hash.
        enriched = json.loads(clean_path.read_text())
        enriched.update(run_id='ABC12345', log_time='2026-07-08 00:00:00',
                        input_id='CHEMBL2971', input_type='chembl', input_source='sql')
        enriched_path = tmp_path / 'enriched.json'
        enriched_path.write_text(json.dumps(enriched))

        p2 = PipeParams.from_json(str(enriched_path))
        assert p2._get_config_hash() == p._get_config_hash()
        for leaked in ('run_id', 'log_time', 'input_id', 'input_type', 'input_source'):
            assert not hasattr(p2, leaked), f"run metadata leaked into params: {leaked}"


# --------------------------------------------------------------------------- #
# Version / hash coupling
# --------------------------------------------------------------------------- #

class TestVersionInHash:

    def test_version_is_single_source(self):
        from molforge import _version
        assert molforge.__version__ == _version.__version__

    def test_hash_depends_on_version(self, monkeypatch):
        p = make_params()
        baseline = p._get_config_hash()
        monkeypatch.setattr(molforge._version, '__version__', '99.9.9')
        assert p._get_config_hash() != baseline, "hash must change when the MolForge version changes"

    def test_hash_stable_for_identical_config(self):
        assert make_params()._get_config_hash() == make_params()._get_config_hash()


# --------------------------------------------------------------------------- #
# Derived parameters (re-derived on load, not injected)
# --------------------------------------------------------------------------- #

class TestFlatCurateConfig:
    """Curate serializes flat inputs only; no derived state reaches the config."""

    def _non_default_curate(self):
        return ForgeParams(
            steps=['curate'],
            curate_params={
                'mol_steps': ['removeIsotope', 'removeHs', 'neutralize',
                              'sanitize', 'handleStereo', 'computeProps'],
                'stereo_policy': 'enumerate',
                'assign_policy': 'lowest',
                'max_isomers': 8,
                'random_seed': 7,
            },
        )

    def test_written_block_holds_inputs_only(self, tmp_path):
        p = self._non_default_curate()
        path = tmp_path / 'config.json'
        p.to_json(str(path))
        block = json.loads(path.read_text())['curate']

        assert block['stereo_policy'] == 'enumerate'
        assert block['max_isomers'] == 8
        assert 'stereo_params' not in block
        for flag in ('desalt', 'removeHs', 'canonical', 'kekulize'):
            assert flag not in block

    @pytest.mark.parametrize('ext,writer', [('json', 'to_json'), ('yaml', 'to_yaml')])
    def test_non_default_curate_roundtrips(self, tmp_path, ext, writer):
        p = self._non_default_curate()
        path = tmp_path / f'config.{ext}'
        getattr(p, writer)(str(path))

        reloaded = PipeParams.from_json(str(path)) if ext == 'json' else PipeParams.from_yaml(str(path))
        curate = reloaded.curate_params

        assert curate.stereo_policy == 'enumerate' and curate.max_isomers == 8
        assert curate.mol_steps == p.curate_params.mol_steps
        assert reloaded._get_config_hash() == p._get_config_hash()

    def test_block_level_typo_raises(self):
        with pytest.raises(TypeError, match="Unknown parameter"):
            ForgeParams.from_config({'steps': ['curate'], 'curate': {'mol_stpes': ['sanitize']}})


class TestNestedParams:
    """A params dataclass nested inside another rebuilds recursively."""

    def test_nested_dataclass_field_rebuilt_and_rederived(self):
        outer = _OuterParams(label='y', inner=_InnerParams(scale=3.0))

        # The private derived attribute stays out of the serialized block.
        data = ConfigCodec._to_plain(outer)
        assert data['inner']['scale'] == 3.0
        assert '_scaled' not in data['inner']

        rebuilt = ConfigCodec.build_param_object(_OuterParams, data)
        assert isinstance(rebuilt.inner, _InnerParams), "nested block stayed a dict"
        assert rebuilt.inner.scale == 3.0
        assert rebuilt.inner._scaled == 30, "derived attribute not re-derived"

    def test_pep604_union_nested_field_rebuilt(self):
        # A nested-params field annotated ``X | None`` is still rebuilt recursively.
        outer = _OuterUnionParams(inner=_InnerParams(scale=2.0))
        data = ConfigCodec._to_plain(outer)
        rebuilt = ConfigCodec.build_param_object(_OuterUnionParams, data)
        assert isinstance(rebuilt.inner, _InnerParams), "X | None nested block stayed a dict"
        assert rebuilt.inner._scaled == 20


# --------------------------------------------------------------------------- #
# Chained round-trip across formats
# --------------------------------------------------------------------------- #

class TestChainedRoundTrip:

    def test_hash_stable_across_json_yaml_json_chain(self, tmp_path):
        p0 = make_params()
        h = p0._get_config_hash()

        j1 = tmp_path / 'a.json'
        p0.to_json(str(j1))
        p1 = PipeParams.from_json(str(j1))
        assert p1._get_config_hash() == h

        y1 = tmp_path / 'b.yaml'
        p1.to_yaml(str(y1))
        p2 = PipeParams.from_yaml(str(y1))
        assert p2._get_config_hash() == h

        j2 = tmp_path / 'c.json'
        p2.to_json(str(j2))
        p3 = PipeParams.from_json(str(j2))
        assert p3._get_config_hash() == h

        # The serialized payload is fixed once converted, so the last two files match.
        assert json.loads(j1.read_text()) == json.loads(j2.read_text())


# --------------------------------------------------------------------------- #
# Invalid input handling
# --------------------------------------------------------------------------- #

class TestInvalidInput:

    @pytest.mark.parametrize('ext,content', [
        ('json', '{}'), ('json', 'null'), ('yaml', ''), ('yaml', 'null'),
    ])
    def test_empty_config_file_is_rejected(self, tmp_path, ext, content):
        path = tmp_path / f'empty.{ext}'
        path.write_text(content)
        with pytest.raises(ValueError, match='empty'):
            MolForge(config_path=str(path))

    @pytest.mark.parametrize('ext,content', [('json', '[1, 2]'), ('yaml', '- 1\n- 2')])
    def test_non_mapping_config_is_rejected(self, tmp_path, ext, content):
        path = tmp_path / f'seq.{ext}'
        path.write_text(content)
        with pytest.raises(ValueError, match='mapping'):
            MolForge(config_path=str(path))

    def test_empty_dict_is_valid_programmatically(self):
        # An empty dict is a valid programmatic default; only an empty FILE is
        # rejected. from_config stays lenient so ForgeParams()/build_params() with
        # no config still yield the default pipeline.
        assert ForgeParams.from_config({}).steps == ForgeParams().steps


# --------------------------------------------------------------------------- #
# Unknown keys
# --------------------------------------------------------------------------- #

class TestUnknownKeys:

    def test_unknown_key_warns_and_is_ignored(self):
        # A typo of a step name is dropped (not set) with a warning that suggests
        # the intended key.
        config = {'steps': ['curate'], 'souce': {'backend': 'sql'}}
        with pytest.warns(UserWarning, match="did you mean 'source'"):
            p = ForgeParams.from_config(config)
        assert p.steps == ['curate']

    def test_recognized_config_does_not_warn(self, recwarn):
        ForgeParams.from_config(make_params().to_config())
        assert not any('unrecognized' in str(w.message).lower() for w in recwarn.list)
