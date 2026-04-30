from __future__ import annotations


class HydroMapError(RuntimeError):
    """Base workflow error."""


class ConfigError(HydroMapError):
    """Raised for invalid configuration contracts."""


class StageExecutionError(HydroMapError):
    """Raised when a stage command fails."""


class GPUNotAvailableError(StageExecutionError):
    """Raised when GPU MD is required but CUDA is unavailable."""


class ParityError(HydroMapError):
    """Raised for parity baseline/check failures."""
