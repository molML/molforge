"""Tests for the shared introspection layer."""

import pytest

from molforge.configuration import introspect
from molforge.configuration.params import ForgeParams
from molforge.actors.params import GenerateConfsParams, ChEMBLCuratorParams


class TestAttributeDocs:
    """Test PEP 257 attribute docstring extraction."""

    def test_picks_up_field_docstring(self):
        docs = introspect.attribute_docs(GenerateConfsParams)
        assert 'backend' in docs
        assert 'rdkit' in docs['backend']
        assert 'max_confs' in docs
        assert docs['max_confs']  # non-empty

    def test_inherited_verbose_doc(self):
        # verbose is defined on BaseParams; inheritance should surface it.
        docs = introspect.attribute_docs(GenerateConfsParams)
        assert 'verbose' in docs
        assert docs['verbose']

    def test_robust_on_builtin(self):
        # Source unavailable for builtins -> empty dict, no crash.
        assert introspect.attribute_docs(int) == {}


class TestDescribeParamClass:
    """Test structured field metadata."""

    def test_skips_private_fields(self):
        names = [f['name'] for f in introspect.describe_param_class(ChEMBLCuratorParams)]
        assert not any(n.startswith('_') for n in names)
        assert 'standard_type' in names

    def test_literal_choices_and_bool(self):
        fields = {f['name']: f for f in introspect.describe_param_class(GenerateConfsParams)}

        backend = fields['backend']
        assert backend['choices'] == ['rdkit', 'openeye']
        assert backend['is_bool'] is False
        assert backend['doc']

        dropna = fields['dropna']
        assert dropna['is_bool'] is True
        assert dropna['default'] is True

        max_confs = fields['max_confs']
        assert max_confs['type'] == 'int'
        assert max_confs['default'] == 200

    def test_optional_detection(self):
        fields = {f['name']: f for f in introspect.describe_param_class(ChEMBLCuratorParams)}
        assert fields['standard_units']['is_optional'] is True


class TestDescribeStep:
    """Test full step description."""

    def test_confs_fields_and_backends(self):
        info = introspect.describe_step('confs')
        assert info['step'] == 'confs'
        assert info['is_plugin'] is False
        assert info['param_class'] is GenerateConfsParams

        fields = {f['name']: f for f in info['fields']}
        assert 'backend' in fields
        assert fields['backend']['choices'] == ['rdkit', 'openeye']
        assert fields['backend']['doc']

        # confs registers rdkit / openeye backends
        assert set(info['backends']) == {'rdkit', 'openeye'}

    def test_unknown_step_raises(self):
        with pytest.raises(ValueError):
            introspect.describe_step('not_a_real_step')


class TestListSteps:
    """Test step listing."""

    def test_lists_core_and_plugin(self):
        steps = introspect.list_steps()
        assert 'confs' in steps['core']
        assert 'source' in steps['core']
        assert isinstance(steps['plugin'], list)


class TestForgeParamsIntrospectable:
    """The top-level ForgeParams should also be introspectable."""

    def test_toplevel_fields(self):
        names = [f['name'] for f in introspect.describe_param_class(ForgeParams)]
        assert 'steps' in names
        assert 'output_root' in names


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
