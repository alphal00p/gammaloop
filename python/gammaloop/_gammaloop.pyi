"""
Gammaloop Python API.
"""


from typing import Optional


def cli_wrapper() -> None:
    """ Starts a CLI interface for GammaLoop exposing rust functionalities. """


class NumeratorAwareGroupingOption:

    @classmethod
    def __new__(_cls,
                numerator_aware_grouping_option: Optional[str] = "group_identical_graphs_up_to_scalar_rescaling",
                compare_canonized_numerator: Optional[bool] = True,
                number_of_samples_for_numerator_comparisons: Optional[int] = 5,
                consider_internal_masses_only_in_numerator_isomorphisms: Optional[bool] = True,
                substitute_numerical_masses_when_comparing_numerators: Optional[bool] = False,
                fully_numerical_substitution_when_comparing_numerators: Optional[bool] = True,
                numerical_samples_seed: Optional[int] = 3,
                ) -> NumeratorAwareGroupingOption:
        """ Creates options for grouping diagrams using isomorphisms and taking into account numerator.
        Possible options are: no_grouping, only_detect_zeroes, group_identical_graphs_up_to_sign and group_identical_graphs_up_to_scalar_rescaling
        """


class SnailFilterOptions:

    @classmethod
    def __new__(_cls,
                veto_snails_attached_to_massive_lines: Optional[bool] = False,
                veto_snails_attached_to_massless_lines: Optional[bool] = True,
                veto_only_scaleless_snails: Optional[bool] = False,
                ) -> SnailFilterOptions:
        """ Creates options for vetoing snail diagrams. """


class SewedFilterOptions:
    @classmethod
    def __new__(_cls,
                filter_tadpoles: Optional[bool] = True,
                ) -> SewedFilterOptions:
        """ Creates options for vetoing tadpoles when sewing diagrams. """


class SelfEnergyFilterOptions:
    @classmethod
    def __new__(_cls,
                veto_self_energy_of_massive_lines: Optional[bool] = True,
                veto_self_energy_of_massless_lines: Optional[bool] = True,
                veto_only_scaleless_self_energy: Optional[bool] = False,
                ) -> SelfEnergyFilterOptions:
        """ Creates options for vetoing self-energy diagrams. """


class TadpolesFilterOptions:
    @classmethod
    def __new__(_cls,
                veto_tadpoles_attached_to_massive_lines: Optional[bool] = True,
                veto_tadpoles_attached_to_massless_lines: Optional[bool] = True,
                veto_only_scaleless_tadpoles: Optional[bool] = False,
                veto_cross_section_sewed_tadpoles: Optional[bool] = False,
                ) -> TadpolesFilterOptions:
        """ Creates options for vetoing tadpole diagrams. """


class FeynGenFilters:

    @classmethod
    def __new__(_cls,
                include_external_self_energy: Optional[bool] = False,
                particle_veto: Optional[list[int]] = [],
                max_number_of_bridges: Optional[int] = 0,
                sewed_filter: Optional[SewedFilterOptions] = None,
                self_energy_filter: Optional[SelfEnergyFilterOptions] = None,
                tadpoles_filter: Optional[TadpolesFilterOptions] = None,
                zero_snails_filter: Optional[SnailFilterOptions] = None,
                perturbative_orders: Optional[dict[str, int]] = None,
                coupling_orders: Optional[dict[str,
                                               tuple[int, int | None]]] = {},
                loop_count_range: Optional[tuple[int, int]] = None,
                fermion_loop_count_range: Optional[tuple[int, int]] = None,
                factorized_loop_topologies_count_range: Optional[tuple[int, int]] = None,
                ) -> FeynGenFilters:
        """ Creates a new set of diagram generation filters. """


class FeynGenOptions:

    @classmethod
    def __new__(_cls,
                generation_type: str,
                initial_particles: list[int],
                final_particles: list[list[int]],
                loop_count_range: tuple[int, int],
                cut_blob_range: tuple[int, int],
                cut_spectator_range: tuple[int, int],
                symmetrize_initial_states: bool,
                symmetrize_final_states: bool,
                symmetrize_left_right_states: bool,
                allow_symmetrization_of_external_fermions_in_amplitudes: bool,
                max_multiplicity_for_fast_cut_filter: int,
                amplitude_filters: Optional[FeynGenFilters] = None,
                cross_section_filters: Optional[FeynGenFilters] = None,
                ) -> FeynGenOptions:
        """ Creates options for steering diagram generation.  """


def setup_rust_logging(level: str, format: str) -> None:
    """ Setup logging for the rust backend. """


def atom_to_canonical_string(atom_str: str) -> str:
    """ Converts the string representation of a Symbolica expression to a canonical form. """


class Worker:
    @classmethod
    def __new__(_cls) -> Worker:
        """ Creates a new worker. """

    def load_model(self, file_path: str) -> None:
        """ Loads a model given the file path of its yaml representation. """

    def load_model_from_yaml_str(self, yaml_str: str) -> None:
        """ Loads a model given its yaml string representation. """

    def get_model(self) -> str:
        """ Returns the yaml string representation of the model currently active. """

    def generate_diagrams(self, generation_options: FeynGenOptions,
                          numerator_aware_isomorphism_grouping: NumeratorAwareGroupingOption,
                          filter_self_loop: Optional[bool] = False,
                          graph_prefix: Optional[str] = "GL",
                          selected_graphs: Optional[list[str]] = None,
                          vetoed_graphs: Optional[list[str]] = None,
                          loop_momentum_bases: Optional[dict[str,
                                                             list[str]]] = None,
                          global_prefactor_color: Optional[str] = None,
                          global_prefactor_colorless: Optional[str] = None,
                          num_threads: Optional[int] = None,
                          ) -> list[str]:
        """ Generates diagrams according to the options given in argument and returns their yaml string representation. """

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

    def export_amplitudes(self, export_root: str, amplitude_names: list[str], export_yaml_str: str, no_evaluators: bool) -> None:
        """ Exports the amplitudes given in argument to the export root given in argument, parse export settings as yaml str"""

    def export_expressions(self, export_root: str, amplitude_list: list[str], format: str, export_yaml_str: str) -> None:
        """Exports the numerator and denominator to the export root given in argument in the format which can be 'default' or 'mathematica' or 'latex'."""

    def export_coupling_replacement_rules(self, export_root: str, format: str) -> None:
        """Exports the coupling replacement rules to the export root given in argument. The format can be 'default' or 'mathematica' or 'latex'."""

    def load_amplitude_integrands(self, path_to_settings: str) -> None:
        """ Loads the gammalooop integrand objects"""

    def write_default_settings(self, export_root: str) -> None:
        """ Writes the default settings to the export root given in argument. """

    def inspect_integrand(self, integrand: str, pt: list[float], term: list[int], force_radius: bool, is_momentum_space: bool, use_f128: bool) -> tuple[float, float]:
        """ Inspects the integrand given in argument at the point given in argument. """

    def inspect_lmw_integrand(self, integrand: str, workspace_path: str, use_f128: bool) -> tuple[float, float]:
        """ inspects the integrand given at the max weight point of the previous run. """

    def integrate_integrand(self, integrand: str, num_cores: int, result_path: str, workspace_path: str, target: tuple[float, float] | None) -> list[tuple[float, float]]:
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

    def sync(self) -> None:
        """sync the worker"""
