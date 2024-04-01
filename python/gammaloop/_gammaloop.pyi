"""
Gammaloop Python API.
"""


def cli_wrapper() -> None:
    """ Starts a CLI interface for GammaLoop exposing rust functionalities. """


class Worker:
    @classmethod
    def new(_cls) -> Worker:
        """ Creates a new worker. """

    def load_model(self, file_path: str) -> None:
        """ Loads a model given the file path of its yaml representation. """

    def load_model_from_yaml_str(self, yaml_str: str) -> None:
        """ Loads a model given its yaml string representation. """

    def get_model(self) -> str:
        """ Returns the yaml string representation of the model currently active. """

    def add_cross_section_from_yaml_str(self, yaml_str: str) -> None:
        """ Adds a cross section to the internal list of the worker given its yaml string representation. """

    def load_cross_sections(self, file_path: str) -> None:
        """ Loads a list of cross sections from a yaml file specifying them located at the path given in argument. """

    def load_cross_sections_from_yaml_str(self, yaml_str: str) -> None:
        """ Loads a list of cross sections from a yaml string representation of that list. """

    def get_cross_sections(self) -> str:
        """ Returns the yaml string representation of the list of all cross sections currently loaded in the worker. """

    def reset_cross_sections(self) -> None:
        """ Resets the internal list of cross sections of the worker. """

    def add_amplitude_from_yaml_str(self, yaml_str: str) -> None:
        """ Adds an amplitude to the internal list of the worker given its yaml string representation. """

    def load_amplitudes(self, file_path: str) -> None:
        """ Loads a list of amplitudes from a yaml file specifying them located at the path given in argument. """

    def load_amplitudes_from_yaml_str(self, yaml_str: str) -> None:
        """ Loads a list of amplitudes from a yaml string representation of that list. """

    def load_amplitudes_derived_data(self, file_path: str) -> None:
        """ Loads the derived data associated to the amplitudes. """

    def get_amplitudes(self) -> str:
        """ Returns the yaml string representation of the list of all amplitudes currently loaded in the worker. """

    def reset_amplitudes(self) -> None:
        """ Resets the internal list of amplitudes of the worker. """

    def export_cross_sections(self, export_root: str, cross_section_names: list[str]) -> None:
        """ Exports the cross sections given in argument to the export root given in argument. """

    def export_amplitudes(self, export_root: str, amplitude_names: list[str]) -> None:
        """ Exports the amplitudes given in argument to the export root given in argument. """

    def export_expressions(self,export_root:str,format:str)->None:
        """Exports the numerator and denominator to the export root given in argument in the format which can be 'default' or 'mathematica' or 'latex'."""


    def export_coupling_replacement_rules(self,export_root:str,format:str)->None:
        """Exports the coupling replacement rules to the export root given in argument. The format can be 'default' or 'mathematica' or 'latex'."""

    def load_amplitude_integrands(self, path_to_settings: str) -> None:
        """ Loads the gammalooop integrand objects"""

    def write_default_settings(self, export_root: str) -> None:
        """ Writes the default settings to the export root given in argument. """

    def inspect_integrand(self, integrand: str, pt: list[float], term: list[int], force_radius: bool, is_momentum_space: bool, use_f128: bool) -> None:
        """ Inspects the integrand given in argument at the point given in argument. """

    def integrate_integrand(self, integrand: str, num_cores: int, result_path: str, workspace_path: str, target: tuple[float, float] | None) -> None:
        """ Integrates the integrand given in argument over the target given in argument. """

    def load_master_node(self, integrand: str) -> None:
        """Setup a master node for the integrand given in argument."""

    def write_batch_input(self, num_cores: int, num_samples: int, export_grid: bool, output_accumulator: bool, workspace_path: str, job_id: int) -> None:
        """ Writes a batch input file for a integration job """

    def process_batch_output(self, workspace_path: str, job_id: int) -> None:
        """process the output of a job"""

    def display_master_node_status(self) -> None:
        """display information about the current run"""

    def update_iter(self) -> None:
        """finish the iteration"""
