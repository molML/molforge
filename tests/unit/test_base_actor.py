# Tests for new BaseActor class

import pytest
import pandas as pd
from molforge.actors.base import BaseActor
from molforge.actors.protocol import ActorInput, ActorOutput
from molforge.configuration.context import PipelineContext
from molforge.actors.params.base import BaseParams
from dataclasses import dataclass


# Create a simple test actor for testing
@dataclass
class TestActorParams(BaseParams):
    test_param: str = "default"

    def _validate_params(self):
        pass


class TestActor(BaseActor):
    """Minimal test actor implementation."""

    @property
    def _nametag(self) -> str:
        return "[TEST]"

    @property
    def required_columns(self):
        return ['input_col']

    @property
    def output_columns(self):
        return ['output_col']

    def process(self, data: pd.DataFrame) -> pd.DataFrame:
        """Simple processing: add output column."""
        data = data.copy()
        data['output_col'] = data['input_col'] * 2
        return data


@pytest.fixture
def output_dir(tmp_path):
    """Unified output directory for all test artifacts."""
    out_dir = tmp_path / "test_output"
    out_dir.mkdir(exist_ok=True)
    return out_dir

def test_base_actor_init(output_dir):
    """Test BaseActor initialization."""
    params = TestActorParams(test_param="value")

    actor = TestActor(params, logger=None)

    assert actor._params is params
    assert actor._context is None  # Context set during execution
    assert actor.test_param == "value"  # From _init_from_params


def test_base_actor_validate_input():
    """Test input validation."""
    params = TestActorParams()
    actor = TestActor(params)

    # Valid input
    df_valid = pd.DataFrame({'input_col': [1, 2, 3]})
    assert actor.validate_input(df_valid) is True

    # Invalid input (missing column)
    df_invalid = pd.DataFrame({'wrong_col': [1, 2, 3]})
    assert actor.validate_input(df_invalid) is False



def test_base_actor_call_success(output_dir):
    """Test successful actor execution."""
    params = TestActorParams(verbose=False)
    context = PipelineContext(run_id='test', output_dir=str(output_dir))
    actor = TestActor(params, context)

    df = pd.DataFrame({'input_col': [1, 2, 3]})
    actor_input = ActorInput(data=df, context=context)

    result = actor(actor_input)

    assert isinstance(result, ActorOutput)
    assert result.success is True
    assert 'output_col' in result.data.columns
    assert list(result.data['output_col']) == [2, 4, 6]



def test_base_actor_call_validation_failure(output_dir):
    """Test actor execution with validation failure."""
    params = TestActorParams(verbose=False)
    context = PipelineContext(run_id='test', output_dir=str(output_dir))
    actor = TestActor(params, logger=None)

    # Missing required column
    df = pd.DataFrame({'wrong_col': [1, 2, 3]})
    actor_input = ActorInput(data=df, context=context)

    result = actor(actor_input)

    assert isinstance(result, ActorOutput)
    assert result.success is False
    assert 'error' in result.metadata



def test_base_actor_context_access(output_dir):
    """Test context get/set methods."""
    params = TestActorParams()
    context = PipelineContext(run_id='test', output_dir=str(output_dir))
    actor = TestActor(params, logger=None)

    # Set context by calling actor (which sets _context)
    df = pd.DataFrame({'input_col': [1, 2, 3]})
    actor_input = ActorInput(data=df, context=context)
    actor(actor_input)

    # Set context value
    actor.set_context('test_key', 'test_value')

    # Get it back
    assert actor.get_context('test_key') == 'test_value'

    # Get with default
    assert actor.get_context('missing', 'default') == 'default'
