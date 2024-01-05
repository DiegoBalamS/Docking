from importlib.metadata import PackageNotFoundError, version

from ligprepare import prep

from target import prept

from .target import list_all_target_names, load_target
from .utils import (
    'load_target', 'list_all_target_names', 'setup_logger', 'DockstringError', 'CanonicalizationError', 'ParsingError',
    'OutputError', 'SanityError', 'EmbeddingError', 'StructureOptimizationError', 'FormatConversionError',
    'ProtonationError', 'PoseProcessingError', 'VinaError', 'DockingError'
)

try:
    __version__ = version("dockstring")
except PackageNotFoundError:
    # package is not installed
    pass
