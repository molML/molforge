"""Tests for the interactive configuration wizard (non-interactive, mocked)."""

from unittest.mock import MagicMock

from molforge.cli import wizard


def _prompt(value):
    """A fake questionary prompt whose ``.unsafe_ask()`` returns ``value``."""
    prompt = MagicMock()
    prompt.unsafe_ask.return_value = value
    return prompt


class TestDefaultSummary:
    """The per-step default-summary helper."""

    def test_renders_name_default_pairs(self):
        fields_meta = [
            {"name": "backend", "default": "rdkit"},
            {"name": "max_confs", "default": 200},
            {"name": "rms_threshold", "default": 0.5},
        ]
        summary = wizard._default_summary(fields_meta)
        assert summary == "backend=rdkit, max_confs=200, rms_threshold=0.5"

    def test_empty_fields(self):
        assert wizard._default_summary([]) == ""


class TestFieldOverride:
    """The field -> override string helper."""

    def test_choice_changed_emits_override(self):
        q = MagicMock()
        q.select.return_value = _prompt("api")
        meta = {"name": "backend", "default": "sql", "choices": ["sql", "api"]}
        assert wizard._field_override(q, "source.", meta) == "source.backend=api"

    def test_choice_unchanged_no_override(self):
        q = MagicMock()
        q.select.return_value = _prompt("sql")
        meta = {"name": "backend", "default": "sql", "choices": ["sql", "api"]}
        assert wizard._field_override(q, "source.", meta) is None

    def test_bool_changed_emits_override(self):
        q = MagicMock()
        q.confirm.return_value = _prompt(False)
        meta = {"name": "dropna", "default": True, "is_bool": True, "choices": None}
        assert wizard._field_override(q, "curate.", meta) == "curate.dropna=false"

    def test_bool_unchanged_no_override(self):
        q = MagicMock()
        q.confirm.return_value = _prompt(True)
        meta = {"name": "dropna", "default": True, "is_bool": True, "choices": None}
        assert wizard._field_override(q, "curate.", meta) is None

    def test_text_changed_emits_override(self):
        q = MagicMock()
        q.text.return_value = _prompt("500")
        meta = {"name": "max_confs", "default": 200, "type": "int", "choices": None}
        assert wizard._field_override(q, "confs.", meta) == "confs.max_confs=500"

    def test_text_unchanged_no_override(self):
        q = MagicMock()
        q.text.return_value = _prompt("200")
        meta = {"name": "max_confs", "default": 200, "type": "int", "choices": None}
        assert wizard._field_override(q, "confs.", meta) is None

    def test_toplevel_prefix_empty(self):
        q = MagicMock()
        q.text.return_value = _prompt("./results")
        meta = {"name": "output_root", "default": "MOLFORGE_OUT", "type": "str", "choices": None}
        assert wizard._field_override(q, "", meta) == "output_root=./results"

    def test_complex_blank_no_override(self):
        q = MagicMock()
        q.text.return_value = _prompt("")
        meta = {"name": "mapping", "default": {}, "type": "Dict", "choices": None}
        assert wizard._field_override(q, "step.", meta) is None

    def test_complex_json_emits_override(self):
        q = MagicMock()
        q.text.return_value = _prompt('{"a": 1}')
        meta = {"name": "mapping", "default": {}, "type": "Dict", "choices": None}
        assert wizard._field_override(q, "step.", meta) == 'step.mapping={"a": 1}'


class TestRunWizardCancel:
    """Ctrl+C at any prompt cancels the wizard."""

    def test_keyboard_interrupt_returns_none(self, monkeypatch, capsys):
        q = MagicMock()
        q.checkbox.return_value.unsafe_ask.side_effect = KeyboardInterrupt
        monkeypatch.setattr(wizard, "_require_questionary", lambda: q)

        result = wizard.run_wizard()
        assert result is None
        assert "cancelled" in capsys.readouterr().out.lower()


class TestQuestionaryAvailable:
    """The availability probe."""

    def test_returns_bool(self):
        assert isinstance(wizard.questionary_available(), bool)
