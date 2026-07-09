"""
Centralized step name constants for type-safe actor references.

Using constants instead of string literals provides:
- IDE autocomplete
- Refactoring safety
- Typo prevention
- Single source of truth
"""

from typing import List


class Steps:
    """
    Step name constants for all core actors.

    These constants are used throughout the codebase for:
    - Pipeline step configuration
    - Actor dependencies
    - Context actor/result retrieval
    - Test assertions

    Example:
        >>> from molforge.configuration.steps import Steps
        >>> params = ForgeParams(steps=[Steps.SOURCE, Steps.CHEMBL, Steps.CURATE])
        >>> gc_actor = context.get_actor(Steps.CONFS)
    """

    # Data Source
    SOURCE = 'source'

    # Curation & Processing
    CHEMBL = 'chembl'
    CURATE = 'curate'
    TOKENS = 'tokens'
    DISTRIBUTIONS = 'distributions'

    # 3D Structure
    CONFS = 'confs'

    @classmethod
    def all(cls) -> List[str]:
        """
        Get list of all core step names.

        Returns:
            List of step name strings
        """
        return [
            value for name, value in vars(cls).items()
            if not name.startswith('_') and isinstance(value, str) and name.isupper()
        ]

    @classmethod
    def validate(cls, step_name: str) -> bool:
        """
        Check if a step name is a valid core step.

        Args:
            step_name: Step name to validate

        Returns:
            True if step_name is in core steps
        """
        return step_name in cls.all()

    @classmethod
    def max_length(cls) -> int:
        """Get the length of the longest step name."""
        return max(len(s) for s in cls.all())