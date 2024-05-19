DEFAULT_EXPTIME = 30.0

# use lower Y to simplify code
FOCUS_DEPTH_DICT = {'u': 0, 'g': 0, 'r': 0, 'i': 0, 'z': 0, 'y': -0.6}

PIXEL_SCALE = 0.2
DEFAULT_COSMIC_RAY_RATE = 0.2
MIN_STAMP_SIZE = 32
MAX_STAMP_SIZE = 4096

# fluxes higher than this are draw with FFT
# see get_initial_draw_method
FFT_FLUX_THRESH = 1.e6
FFT_SB_THRESH = 2.0e5

DEFAULT_MAX_FLUX_SIMPLE = 100
DEFAULT_MINMAG = -1000
DEFAULT_MAXN = 1_000_000
