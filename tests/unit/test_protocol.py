# Tests for new actor protocol components

import pytest
import pandas as pd
from molforge.actors.protocol import ActorInput, ActorOutput
from molforge.configuration.context import PipelineContext


def test_actor_input_creation():
    """Test ActorInput creation."""
    df = pd.DataFrame({'col': [1, 2, 3]})
    context = PipelineContext(run_id='test', output_dir='/tmp')

    actor_input = ActorInput(data=df, context=context)

    assert actor_input.data is df
    assert actor_input.context is context


def test_actor_output_creation():
    """Test ActorOutput creation with defaults."""
    df = pd.DataFrame({'result': [1, 2, 3]})

    output = ActorOutput(data=df)

    assert output.data is df
    assert output.success is True
    assert output.metadata == {}
    assert output.endpoint is None


def test_actor_output_with_metadata():
    """Test ActorOutput with custom metadata."""
    df = pd.DataFrame({'result': [1, 2]})
    metadata = {'n_processed': 2, 'source': 'test'}

    output = ActorOutput(data=df, success=True, metadata=metadata)

    assert output.metadata == metadata


def test_pipeline_context_shared_state():
    """Test PipelineContext state management."""
    context = PipelineContext(run_id='test123', output_dir='/tmp')

    # Set and get
    context.set('key1', 'value1')
    assert context.get('key1') == 'value1'

    # Get with default
    assert context.get('missing', 'default') == 'default'


def test_pipeline_context_actor_results():
    """Test PipelineContext actor result storage."""
    context = PipelineContext(run_id='test', output_dir='/tmp')
    df = pd.DataFrame({'col': [1]})
    result = ActorOutput(data=df, success=True)

    # Store actor result
    context.set_actor_result('GC', result)

    # Retrieve it
    retrieved = context.get_actor_result('GC')
    assert retrieved is result
    assert retrieved.success is True
