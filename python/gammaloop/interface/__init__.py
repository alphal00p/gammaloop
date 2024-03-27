from gammaloop import check_gammaloop_dependencies

# Make sure that when gammaloop is imported as a library, the sanity checks are run
# Do not place this check in the top-level __init__.py of the gammaloop module because
# we need to keep it possible for the cli() function to first issue calls to build_dependencies.sh
check_gammaloop_dependencies()
