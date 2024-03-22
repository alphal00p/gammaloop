from gammaloop.tests.common import pjoin
from gammaloop.misc.common import GL_PATH
from gammaloop.base_objects.param_card import ParamCard, ParamCardWriter
from gammaloop.base_objects.model import InputParamCard, Model
from pathlib import Path
import pytest
# from gammaloop.misc.utils import expression_to_string
# from pprint import pprint, pformat


class TestParamCard:

    def test_load_sm_param_card(self, tmpdir_factory: pytest.TempPathFactory):

        SM_PARAMETERS = {'aewm1', 'gf', 'as', 'mc', 'mb', 'mt', 'me', 'mm', 'mta', 'mz', 'wtau', 'wt', 'wh', 'ww', 'wz',
                         'mh', 'lamws', 'aws', 'rhows', 'etaws', 'ymc', 'ymb', 'ymt', 'yme', 'ymm', 'ymtau'}

        for restrict_card, (target_n_parameters, target_n_interactions, target_n_couplings) in [
            ('restrict_default.dat', (15, 119, 72)),
            ('restrict_no_widths.dat', (11, 119, 72)),
            ('restrict_lepton_masses.dat', (20, 127, 80)),
            ('restrict_c_mass.dat', (17, 123, 76)),
            ('restrict_no_b_mass.dat', (13, 117, 68)),
            ('restrict_zeromass_ckm.dat', (15, 129, 80)),
            ('restrict_ckm.dat', (19, 139, 92)),
            ('restrict_no_masses.dat', (11, 113, 64)),
            ('restrict_no_tau_mass.dat', (13, 115, 68)),
            (None, (26, 153, 108))
        ]:
            model = Model.from_ufo_model('sm')

            if restrict_card is not None:
                param_card = ParamCard(
                    pjoin(GL_PATH, 'data', 'models', 'sm', restrict_card))
                input_param_card = InputParamCard.from_param_card(
                    param_card, model)
            else:
                input_param_card = InputParamCard.default_from_model(model)
                model.apply_input_param_card(input_param_card, simplify=False)
                default_card_path = Path(tmpdir_factory.mktemp(
                    "default_param_card_tmp_directory")).joinpath("sm_model_default_lha_param_card.dat")
                ParamCardWriter.write(
                    default_card_path, model, generic=True)
                param_card = ParamCard(str(default_card_path))

            all_params, _restrictions = param_card.analyze_param_card(model)
            assert set(all_params.keys()) == SM_PARAMETERS

            model.apply_input_param_card(input_param_card, simplify=True)
            assert len(input_param_card) == target_n_parameters
            assert len(model.vertex_rules) == target_n_interactions
            assert len(model.couplings) == target_n_couplings
            assert not any(p.value is None for p in model.parameters)
            assert not any(c.value is None for c in model.couplings)

    def test_load_scalars_param_card(self, tmpdir_factory: pytest.TempPathFactory):

        SCALARS_PARAMETERS = {'width_scalar_1', 'as',
                              'width_scalar_2', 'aewm1', 'gf', 'mass_scalar_1', 'lam', 'mass_scalar_2'}

        for restrict_card, (target_n_parameters, target_n_interactions, target_n_couplings) in [
            ('restrict_default.dat', (9, 25, 1)),
            (None, (8, 25, 1))
        ]:
            model = Model.from_ufo_model('scalars')

            if restrict_card is not None:
                param_card = ParamCard(
                    pjoin(GL_PATH, 'data', 'models', 'scalars', restrict_card))
                input_param_card = InputParamCard.from_param_card(
                    param_card, model)
            else:
                input_param_card = InputParamCard.default_from_model(model)
                model.apply_input_param_card(input_param_card, simplify=False)
                default_card_path = Path(tmpdir_factory.mktemp(
                    "default_param_card_tmp_directory")).joinpath("scalars_model_default_lha_param_card.dat")
                ParamCardWriter.write(
                    default_card_path, model, generic=True)
                param_card = ParamCard(str(default_card_path))

            all_params, _restrictions = param_card.analyze_param_card(model)
            if restrict_card != 'restrict_default.dat':
                assert set(all_params.keys()) == SCALARS_PARAMETERS
            else:
                assert set(all_params.keys()
                           ) == SCALARS_PARAMETERS ^ {'scalar_0'}

            model.apply_input_param_card(input_param_card, simplify=False)
            assert len(input_param_card) == target_n_parameters
            assert len(model.vertex_rules) == target_n_interactions
            assert len(model.couplings) == target_n_couplings
            assert not any(p.value is None for p in model.parameters)
            assert not any(c.value is None for c in model.couplings)
