"""
Comprehensive test suite for MolForge configuration templates.

Tests all templates for:
- YAML validity
- Parameter correctness
- Schema compliance
- Actor parameter class compatibility
- Configuration loading
"""

import pytest
from pathlib import Path
from typing import Dict, Any
import yaml

from molforge.cli.commands import get_template, load_config
from molforge.configuration.registry import ActorRegistry
from molforge.configuration.steps import Steps


# Template directory
TEMPLATE_DIR = Path(__file__).parent.parent.parent / "molforge" / "cli" / "templates"

# All available templates
AVAILABLE_TEMPLATES = [
    "basic",
    "api", 
    "conformers",
    "distributions",
    "full"
]

# Expected steps for each template
EXPECTED_STEPS = {
    "basic": ["source", "chembl", "curate"],
    "api": ["source", "chembl", "curate"],
    "conformers": ["source", "chembl", "curate", "confs"],
    "distributions": ["source", "chembl", "curate", "tokens", "distributions"],
    "full": ["source", "chembl", "curate", "tokens", "distributions", "confs"],
}

# Expected backend for each template
EXPECTED_BACKENDS = {
    "basic": "sql",
    "api": "api",
    "conformers": "sql",
    "distributions": "sql",
    "full": "sql",
}


class TestTemplateDiscovery:
    """Test template discovery and loading."""
    
    def test_template_directory_exists(self):
        """Test that template directory exists."""
        assert TEMPLATE_DIR.exists(), f"Template directory not found: {TEMPLATE_DIR}"
        assert TEMPLATE_DIR.is_dir(), f"Template path is not a directory: {TEMPLATE_DIR}"
    
    def test_all_templates_exist(self):
        """Test that all expected templates exist as files."""
        for template_name in AVAILABLE_TEMPLATES:
            template_file = TEMPLATE_DIR / f"{template_name}.yaml"
            assert template_file.exists(), f"Template file not found: {template_file}"
            assert template_file.is_file(), f"Template path is not a file: {template_file}"
    
    def test_get_template_function(self):
        """Test get_template function for all templates."""
        for template_name in AVAILABLE_TEMPLATES:
            template_content = get_template(template_name)
            assert template_content, f"Template {template_name} returned empty content"
            assert isinstance(template_content, str), f"Template {template_name} not a string"
            assert len(template_content) > 0, f"Template {template_name} is empty"
    
    def test_get_template_invalid_name(self):
        """Test get_template with invalid template name."""
        with pytest.raises(ValueError, match="Unknown template"):
            get_template("nonexistent_template")


class TestTemplateYAMLValidity:
    """Test YAML syntax validity of all templates."""
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_yaml_parseable(self, template_name):
        """Test that template YAML can be parsed."""
        template_content = get_template(template_name)
        
        # Should not raise exception
        config = yaml.safe_load(template_content)
        assert config is not None, f"Template {template_name} parsed to None"
        assert isinstance(config, dict), f"Template {template_name} not a dict"
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_yaml_has_version(self, template_name):
        """Test that template has version field."""
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        
        assert "version" in config, f"Template {template_name} missing 'version' field"
        assert config["version"] == "1.0", f"Template {template_name} has wrong version"
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_yaml_has_steps(self, template_name):
        """Test that template has steps field."""
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        
        assert "steps" in config, f"Template {template_name} missing 'steps' field"
        assert isinstance(config["steps"], list), f"Template {template_name} steps not a list"
        assert len(config["steps"]) > 0, f"Template {template_name} has empty steps"
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_yaml_has_output(self, template_name):
        """Test that template has output configuration."""
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        
        assert "output" in config, f"Template {template_name} missing 'output' field"
        assert isinstance(config["output"], dict), f"Template {template_name} output not a dict"
        assert "dir" in config["output"], f"Template {template_name} output missing 'dir'"


class TestTemplateSteps:
    """Test that templates have correct pipeline steps."""
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_expected_steps(self, template_name):
        """Test that template has expected steps."""
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        
        expected = EXPECTED_STEPS[template_name]
        actual = config["steps"]
        
        assert actual == expected, (
            f"Template {template_name} has wrong steps.\n"
            f"Expected: {expected}\n"
            f"Actual:   {actual}"
        )
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_steps_are_valid(self, template_name):
        """Test that all steps in template are valid actor names."""
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        
        # Discover actors
        ActorRegistry.discover_plugins()
        valid_steps = set(Steps.all())
        
        for step in config["steps"]:
            assert step in valid_steps, (
                f"Template {template_name} has invalid step: {step}\n"
                f"Valid steps: {sorted(valid_steps)}"
            )
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_step_configs_exist(self, template_name):
        """Test that each step has a configuration section."""
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        
        for step in config["steps"]:
            assert step in config, (
                f"Template {template_name} missing config section for step: {step}"
            )
            assert isinstance(config[step], dict), (
                f"Template {template_name} config for {step} is not a dict"
            )


class TestTemplateBackends:
    """Test that templates use correct backends."""
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_source_backend(self, template_name):
        """Test that source backend is correctly set."""
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        
        expected_backend = EXPECTED_BACKENDS[template_name]
        
        if "source" in config["steps"]:
            assert "backend" in config["source"], (
                f"Template {template_name} source missing 'backend' field"
            )
            actual_backend = config["source"]["backend"]
            assert actual_backend == expected_backend, (
                f"Template {template_name} has wrong backend.\n"
                f"Expected: {expected_backend}\n"
                f"Actual:   {actual_backend}"
            )


class TestTemplateParameterCompatibility:
    """Test that template parameters match actor parameter classes."""
    
    @pytest.fixture
    def actor_registry(self):
        """Initialize actor registry."""
        ActorRegistry.discover_plugins()
        return ActorRegistry
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_source_params(self, template_name, actor_registry):
        """Test that source parameters are valid."""
        if "source" not in EXPECTED_STEPS[template_name]:
            pytest.skip(f"Template {template_name} doesn't use source")
        
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        source_config = config.get("source", {})
        
        # Get parameter class
        from molforge.actors.params.source import ChEMBLSourceParams
        from dataclasses import fields
        
        valid_params = {f.name for f in fields(ChEMBLSourceParams)}
        
        for param_name in source_config.keys():
            assert param_name in valid_params, (
                f"Template {template_name} has invalid source param: {param_name}\n"
                f"Valid params: {sorted(valid_params)}"
            )
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_chembl_params(self, template_name, actor_registry):
        """Test that chembl parameters are valid."""
        if "chembl" not in EXPECTED_STEPS[template_name]:
            pytest.skip(f"Template {template_name} doesn't use chembl")
        
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        chembl_config = config.get("chembl", {})
        
        # Get parameter class
        from molforge.actors.params.chembl import ChEMBLCuratorParams
        from dataclasses import fields
        
        valid_params = {f.name for f in fields(ChEMBLCuratorParams)}
        
        for param_name in chembl_config.keys():
            assert param_name in valid_params, (
                f"Template {template_name} has invalid chembl param: {param_name}\n"
                f"Valid params: {sorted(valid_params)}"
            )
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_curate_params(self, template_name, actor_registry):
        """Test that curate parameters are valid."""
        if "curate" not in EXPECTED_STEPS[template_name]:
            pytest.skip(f"Template {template_name} doesn't use curate")
        
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        curate_config = config.get("curate", {})
        
        # Get parameter class
        from molforge.actors.params.curate import CurateMolParams
        from dataclasses import fields
        
        valid_params = {f.name for f in fields(CurateMolParams)}
        
        for param_name in curate_config.keys():
            assert param_name in valid_params, (
                f"Template {template_name} has invalid curate param: {param_name}\n"
                f"Valid params: {sorted(valid_params)}"
            )
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_tokens_params(self, template_name, actor_registry):
        """Test that tokens parameters are valid."""
        if "tokens" not in EXPECTED_STEPS[template_name]:
            pytest.skip(f"Template {template_name} doesn't use tokens")
        
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        tokens_config = config.get("tokens", {})
        
        # Get parameter class
        from molforge.actors.params.tokens import TokenizeDataParams
        from dataclasses import fields
        
        valid_params = {f.name for f in fields(TokenizeDataParams)}
        
        for param_name in tokens_config.keys():
            assert param_name in valid_params, (
                f"Template {template_name} has invalid tokens param: {param_name}\n"
                f"Valid params: {sorted(valid_params)}"
            )
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_distributions_params(self, template_name, actor_registry):
        """Test that distributions parameters are valid."""
        if "distributions" not in EXPECTED_STEPS[template_name]:
            pytest.skip(f"Template {template_name} doesn't use distributions")
        
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        dist_config = config.get("distributions", {})
        
        # Get parameter class
        from molforge.actors.params.distributions import CurateDistributionParams
        from dataclasses import fields
        
        valid_params = {f.name for f in fields(CurateDistributionParams)}
        
        for param_name in dist_config.keys():
            assert param_name in valid_params, (
                f"Template {template_name} has invalid distributions param: {param_name}\n"
                f"Valid params: {sorted(valid_params)}"
            )
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_confs_params(self, template_name, actor_registry):
        """Test that confs parameters are valid."""
        if "confs" not in EXPECTED_STEPS[template_name]:
            pytest.skip(f"Template {template_name} doesn't use confs")
        
        template_content = get_template(template_name)
        config = yaml.safe_load(template_content)
        confs_config = config.get("confs", {})
        
        # Get parameter class
        from molforge.actors.params.confs import GenerateConfsParams
        from dataclasses import fields
        
        valid_params = {f.name for f in fields(GenerateConfsParams)}
        
        for param_name in confs_config.keys():
            assert param_name in valid_params, (
                f"Template {template_name} has invalid confs param: {param_name}\n"
                f"Valid params: {sorted(valid_params)}"
            )


class TestTemplateConfigLoading:
    """Test that templates can be loaded as configurations."""
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_template_loads_as_config(self, template_name, tmp_path):
        """Test that template can be loaded as a valid configuration."""
        # Get template content
        template_content = get_template(template_name)
        
        # Write to temporary file
        config_file = tmp_path / f"{template_name}.yaml"
        config_file.write_text(template_content)
        
        # Load configuration
        # This should not raise an exception
        params = load_config(str(config_file))
        
        # Verify basic properties
        assert params is not None, f"Template {template_name} loaded as None"
        assert hasattr(params, 'steps'), f"Template {template_name} missing steps"
        assert params.steps == EXPECTED_STEPS[template_name], (
            f"Template {template_name} loaded with wrong steps"
        )


class TestTemplateSpecificRequirements:
    """Test template-specific requirements and characteristics."""
    
    def test_basic_template_minimal(self):
        """Test that basic template is minimal and simple."""
        config = yaml.safe_load(get_template("basic"))
        
        # Should have exactly 3 steps
        assert len(config["steps"]) == 3
        
        # Should use SQL backend
        assert config["source"]["backend"] == "sql"
        
        # Should have basic curation
        assert "mol_steps" in config["curate"]
        assert len(config["curate"]["mol_steps"]) <= 4  # Keep it simple
    
    def test_api_template_uses_api(self):
        """Test that api template uses API backend."""
        config = yaml.safe_load(get_template("api"))
        
        # Must use API backend
        assert config["source"]["backend"] == "api"
        
        # Should have minimal steps (same as basic)
        assert len(config["steps"]) == 3
    
    def test_conformers_template_has_confs(self):
        """Test that conformers template includes conformer generation."""
        config = yaml.safe_load(get_template("conformers"))
        
        # Must include confs step
        assert "confs" in config["steps"]
        
        # Must have confs configuration
        assert "confs" in config
        assert "backend" in config["confs"]
        assert "max_confs" in config["confs"]
        assert "rms_threshold" in config["confs"]
    
    def test_distributions_template_has_filtering(self):
        """Test that distributions template includes distribution filtering."""
        config = yaml.safe_load(get_template("distributions"))
        
        # Must include distributions step
        assert "distributions" in config["steps"]
        
        # Must have distributions configuration
        assert "distributions" in config
        
        # Should have some threshold configuration
        dist_config = config["distributions"]
        has_threshold = (
            "global_statistical_threshold" in dist_config or
            "global_quantile_threshold" in dist_config or
            "thresholds" in dist_config
        )
        assert has_threshold, "distributions template missing threshold configuration"
    
    def test_full_template_has_all_steps(self):
        """Test that full template includes all major steps."""
        config = yaml.safe_load(get_template("full"))
        
        # Should have most steps (at least 5)
        assert len(config["steps"]) >= 5
        
        # Must include key steps
        required_steps = ["source", "chembl", "curate", "tokens", "distributions", "confs"]
        for step in required_steps:
            assert step in config["steps"], f"full template missing required step: {step}"
        
        # Should have checkpointing enabled
        assert config["output"].get("checkpoints") == True


class TestTemplateConsistency:
    """Test consistency across templates."""
    
    def test_all_use_version_1(self):
        """Test that all templates declare the v1.0 config version."""
        for template_name in AVAILABLE_TEMPLATES:
            config = yaml.safe_load(get_template(template_name))
            assert config["version"] == "1.0", (
                f"Template {template_name} not using version 1.0"
            )
    
    def test_all_have_output_dir(self):
        """Test that all templates specify output directory."""
        for template_name in AVAILABLE_TEMPLATES:
            config = yaml.safe_load(get_template(template_name))
            assert "dir" in config["output"], (
                f"Template {template_name} missing output.dir"
            )
            assert isinstance(config["output"]["dir"], str), (
                f"Template {template_name} output.dir not a string"
            )
    
    def test_sql_templates_have_auto_download(self):
        """Test that SQL templates have auto_download configured."""
        for template_name in AVAILABLE_TEMPLATES:
            config = yaml.safe_load(get_template(template_name))
            
            if config.get("source", {}).get("backend") == "sql":
                assert "auto_download" in config["source"], (
                    f"SQL template {template_name} missing auto_download"
                )
    
    def test_curate_has_stereo_policy(self):
        """Test that curate step has stereochemistry policy."""
        for template_name in AVAILABLE_TEMPLATES:
            config = yaml.safe_load(get_template(template_name))
            
            if "curate" in config["steps"]:
                assert "stereo_policy" in config["curate"], (
                    f"Template {template_name} curate missing stereo_policy"
                )


class TestTemplateDocumentation:
    """Test that templates have proper documentation/comments."""
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_template_has_header_comment(self, template_name):
        """Test that template has header comment."""
        template_content = get_template(template_name)
        
        # Should start with comment
        lines = template_content.strip().split('\n')
        assert lines[0].startswith('#'), (
            f"Template {template_name} missing header comment"
        )
    
    @pytest.mark.parametrize("template_name", AVAILABLE_TEMPLATES)
    def test_template_has_section_comments(self, template_name):
        """Test that template has section comments."""
        template_content = get_template(template_name)
        
        # Should have at least 3 comment lines
        comment_lines = [line for line in template_content.split('\n') if line.strip().startswith('#')]
        assert len(comment_lines) >= 3, (
            f"Template {template_name} has too few comments (less than 3)"
        )


# Pytest marks for categorization
pytestmark = [
    pytest.mark.cli,
    pytest.mark.templates,
]


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "--tb=short"])
