# Planned features and interactions

All configuration settings (not process specific) should be settable from a file (a default is provided). It should be overridable at run time. Overriding does not modify the file.

The model file uses the ufo format. When loading, this gets converted to a yaml file. All _independent_ values (external parameters as per UFO) in this file should be overridable, and all _dependent_ values are computed and modified when the values are updated. This modifies the yaml file. When using any parameter, or derived value we avoid recomputing it and just use its value.