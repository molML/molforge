"""
Tests for the shared parameter builder.

None of these tests may trigger a ChEMBL download: pipelines either omit the
'source' step or force the API backend (which performs no local DB setup) via
overrides / config, keeping the suite free of network I/O.
"""

import json

import pytest

from molforge.cli.builder import build_params, parse_override
from molforge.configuration.params import ForgeParams


class TestBuildWithSteps:
    """Steps are honoured and every selected step's params get initialized."""

    def test_non_default_subset_initializes_confs(self):
        # 'confs' is NOT in the default pipeline; the old CLI mutated .steps
        # after construction so confs_params stayed None. Building with steps up
        # front must initialize AND validate confs_params.
        # (source is forced onto the API backend so no ChEMBL DB download occurs.)
        params = build_params(
            steps=['source', 'curate', 'confs'],
            overrides=['source.backend=api'],
        )
        assert params.steps == ['source', 'curate', 'confs']
        assert params.confs_params is not None
        assert params.confs_params.backend == 'rdkit'
        # A step left out of the pipeline is not initialized.
        assert params.tokens_params is None

    def test_steps_accepts_comma_string(self):
        params = build_params(steps='curate,tokens')
        assert params.steps == ['curate', 'tokens']

    def test_default_when_no_steps(self):
        params = build_params(steps=['curate'])
        assert params.steps == ['curate']


class TestBuildWithConfig:
    """Config loading plus steps/overrides applied on top of it."""

    def _write_config(self, tmp_path):
        config = tmp_path / "config.yaml"
        config.write_text(
            """
steps:
  - source
  - chembl
  - curate

source:
  backend: api
  n: 5000

chembl:
  standard_type: IC50
  standard_units: nM

output:
  dir: ./from_config
  checkpoints: false
"""
        )
        return config

    def test_config_loaded(self, tmp_path):
        params = build_params(config_path=str(self._write_config(tmp_path)))
        assert params.steps == ['source', 'chembl', 'curate']
        assert params.source_params.backend == 'api'
        assert params.source_params.n == 5000
        assert params.output_root == './from_config'

    def test_steps_override_config(self, tmp_path):
        # --steps must win over the config's steps (old CLI ignored it).
        params = build_params(
            config_path=str(self._write_config(tmp_path)),
            steps=['curate', 'tokens'],
        )
        assert params.steps == ['curate', 'tokens']

    def test_overrides_applied_on_top_of_config(self, tmp_path):
        params = build_params(
            config_path=str(self._write_config(tmp_path)),
            overrides=['source.n=42', 'chembl.standard_type=Ki'],
        )
        assert params.source_params.n == 42
        assert params.chembl_params.standard_type == 'Ki'
        # Untouched config values are preserved.
        assert params.source_params.backend == 'api'

    def test_json_config(self, tmp_path):
        config = tmp_path / "config.json"
        config.write_text(json.dumps({
            "steps": ["curate", "tokens"],
            "output": {"dir": "./j"},
        }))
        params = build_params(config_path=str(config))
        assert params.steps == ['curate', 'tokens']
        assert params.output_root == './j'


class TestOverrideCoercion:
    """Type coercion + step routing for --set overrides."""

    def test_int_coercion(self):
        params = build_params(steps=['confs'], overrides=['confs.max_confs=50'])
        assert params.confs_params.max_confs == 50
        assert isinstance(params.confs_params.max_confs, int)

    def test_float_coercion(self):
        params = build_params(steps=['confs'], overrides=['confs.rms_threshold=0.25'])
        assert params.confs_params.rms_threshold == 0.25

    def test_bool_coercion(self):
        params = build_params(steps=['confs'], overrides=['confs.dropna=false'])
        assert params.confs_params.dropna is False

    def test_literal_choice_validation(self):
        params = build_params(steps=['confs'], overrides=['confs.backend=openeye'])
        assert params.confs_params.backend == 'openeye'

    def test_literal_choice_step_routing_chembl(self):
        params = build_params(
            steps=['source', 'chembl', 'curate'],
            overrides=['chembl.standard_type=Ki', 'source.backend=api'],
        )
        assert params.chembl_params.standard_type == 'Ki'
        assert params.source_params.backend == 'api'

    def test_list_coercion(self):
        params = build_params(
            steps=['curate'],
            overrides=['curate.mol_steps=desalt,sanitize'],
        )
        assert params.curate_params.mol_steps == ['desalt', 'sanitize']

    def test_toplevel_override(self):
        params = build_params(steps=['curate'], overrides=['output_root=./custom'])
        assert params.output_root == './custom'

    def test_optional_none_override(self):
        params = build_params(
            steps=['source', 'chembl', 'curate'],
            overrides=['source.backend=api', 'chembl.standard_relation=none'],
        )
        assert params.chembl_params.standard_relation is None


class TestInvalidOverrides:
    """Invalid override handling."""

    def test_missing_equals(self):
        with pytest.raises(ValueError):
            build_params(steps=['confs'], overrides=['confs.max_confs'])

    def test_unknown_step(self):
        with pytest.raises(ValueError):
            build_params(steps=['confs'], overrides=['bogus.field=1'])

    def test_unknown_field(self):
        with pytest.raises(ValueError):
            build_params(steps=['confs'], overrides=['confs.not_a_field=1'])

    def test_bad_int(self):
        with pytest.raises(ValueError):
            build_params(steps=['confs'], overrides=['confs.max_confs=notanint'])

    def test_invalid_literal_choice(self):
        with pytest.raises(ValueError):
            build_params(steps=['confs'], overrides=['confs.backend=amber'])

    def test_unknown_toplevel_key(self):
        with pytest.raises(ValueError):
            build_params(steps=['curate'], overrides=['nonexistent=1'])


class TestParseOverride:
    """Low-level override parsing."""

    def test_toplevel(self):
        assert parse_override('output_root=./x') == (None, 'output_root', './x')

    def test_step_scoped(self):
        assert parse_override('confs.max_confs=50') == ('confs', 'max_confs', '50')

    def test_invalid(self):
        with pytest.raises(ValueError):
            parse_override('nope')


class TestOutputRoot:
    def test_output_root_passthrough(self):
        params = build_params(steps=['curate'], output_root='./out')
        assert params.output_root == './out'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
