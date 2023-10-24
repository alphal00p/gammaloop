#
# Pretty horrenduous code ported from models/check_param_card.py in MadGraph5_aMC@NLO.
# It should be rewritten, but it will do for now
#

from __future__ import annotations
import math
import os
import re
from io import StringIO
from typing import Any
from gammaloop.misc.common import GammaLoopError, logger  # type: ignore
import gammaloop.misc.utils as utils
from pathlib import Path
from typing import TYPE_CHECKING, TextIO
from functools import cmp_to_key
if TYPE_CHECKING:
    from gammaloop.base_objects.model import Model, Particle
    from gammaloop.base_objects.model import Parameter as ModelParameter
pjoin = os.path.join


class InvalidParamCard(GammaLoopError):
    """ An error class for invalid param_card """
    pass


class Parameter(object):
    """A class for a param_card parameter"""

    def __init__(self, block: str = '', lhacode: tuple[int, ...] | None = None, value: float = 0., comment: str = ''):
        """Init the parameter"""

        self.format: str = 'float'

        self.lhablock = block
        if lhacode:
            self.lhacode: tuple[int, ...] = lhacode
        else:
            self.lhacode: tuple[int, ...] = tuple([])
        self.value: float | str = value
        self.comment: str = comment

    def set_block(self, block: str):
        """ set the block name """

        self.lhablock = block

    def load_str(self, text: str):
        """ initialize the information from a str"""

        if '#' in text:
            data, self.comment = text.split('#', 1)
        else:
            data, self.comment = text, ""

        data = data.split()
        if any(d.startswith('scan') for d in data):
            position = [i for i, d in enumerate(
                data) if d.startswith('scan')][0]
            data = data[:position] + [' '.join(data[position:])]
        if not len(data):
            return
        try:
            self.lhacode = tuple([int(d) for d in data[:-1]])
        except Exception:
            self.lhacode = tuple([int(d) for d in data[:-1] if d.isdigit()])
            self.value = ' '.join(data[len(self.lhacode):])
        else:
            self.value = data[-1]

        # convert to number when possible
        try:
            self.value = float(self.value)
        except:
            self.format = 'str'
            pass
        else:
            if self.lhablock == 'modsel':
                self.format = 'int'
                self.value = int(self.value)

    def load_decay(self, text: str):
        """ initialize the decay information from a str"""

        if '#' in text:
            data, self.comment = text.split('#', 1)
        else:
            data, self.comment = text, ""

        if ']]>' in data:
            data = data.split(']]>', 1)[0]

        data = data.split()
        if not len(data):
            return
        lhacode = [int(d) for d in data[2:]]
        lhacode.sort()
        self.lhacode = tuple([len(lhacode)] + lhacode)

        self.value = float(data[0])
        self.format = 'decay_table'

    def __str__(self, precision: str | int = ''):
        """ return a SLAH string """

        format = self.format
        if self.format == 'float':
            try:
                _ = float(self.value)
            except:
                format = 'str'
        self.comment = self.comment.strip()
        if not precision:
            precision = 6

        self.comment = self.comment.strip()
        if format == 'float':
            if self.lhablock == 'decay' and not isinstance(self.value, str):
                return 'DECAY %s %.{0}e # %s'.format(precision) % (' '.join([str(d) for d in self.lhacode]), self.value, self.comment)
            elif self.lhablock == 'decay':
                return 'DECAY %s Auto # %s' % (' '.join([str(d) for d in self.lhacode]), self.comment)
            elif self.lhablock and self.lhablock.startswith('qnumbers'):
                return '      %s %i # %s' % (' '.join([str(d) for d in self.lhacode]), int(self.value), self.comment)
            else:
                return '      %s %.{0}e # %s'.format(precision) % (' '.join([str(d) for d in self.lhacode]), self.value, self.comment)
        elif format == 'int':
            return '      %s %i # %s' % (' '.join([str(d) for d in self.lhacode]), int(self.value), self.comment)
        elif format == 'str':
            if self.lhablock == 'decay':
                return 'DECAY %s %s # %s' % (' '.join([str(d) for d in self.lhacode]), self.value, self.comment)
            return '      %s %s # %s' % (' '.join([str(d) for d in self.lhacode]), self.value, self.comment)
        elif self.format == 'decay_table':
            return '      %e %s # %s' % (self.value, ' '.join([str(d) for d in self.lhacode]), self.comment)
        elif self.format == 'int':
            return '      %s %i # %s' % (' '.join([str(d) for d in self.lhacode]), int(self.value), self.comment)
        else:
            if self.lhablock == 'decay':
                return 'DECAY %s %d # %s' % (' '.join([str(d) for d in self.lhacode]), self.value, self.comment)
            else:
                return '      %s %d # %s' % (' '.join([str(d) for d in self.lhacode]), self.value, self.comment)


class Block(list[Parameter]):
    """ list of parameter """

    def __init__(self, name: str | None = None):
        if name:
            self.name: str = name.lower()
        else:
            self.name: str = ''
        self.scale = None
        self.comment = ''
        self.decay_table: dict[int, Block] = {}
        self.param_dict: dict[tuple[int, ...], Parameter] = {}
        list[Parameter].__init__(self)

    def get(self, lhacode_in: int | tuple[int, ...], default: Any = None) -> Parameter:
        """return the parameter associate to the lhacode"""
        if not self.param_dict:
            self.create_param_dict()

        if isinstance(lhacode_in, int):
            lhacode = (lhacode_in,)
        else:
            lhacode = lhacode_in

        try:
            return self.param_dict[tuple(lhacode)]
        except KeyError:
            if default is None:
                raise KeyError('id %s is not in %s' %
                               (tuple(lhacode), self.name))
            else:
                return Parameter(block=self.name, lhacode=tuple(lhacode), value=default,
                                 comment='not define')

    def remove(self, lhacode: tuple[int, ...]) -> Parameter:  # type: ignore
        """ remove a parameter """
        list[Parameter].remove(self, self.get(lhacode))
        # update the dictionary of key
        return self.param_dict.pop(tuple(lhacode))

    def __eq__(self, other: Block | str, prec: float = 1e-4):  # type: ignore
        """ """

        if isinstance(other, str) and ' ' not in other:
            return self.name.lower() == other.lower()

        if len(self) != len(other):
            return False

        return not any(abs(param.value-other.param_dict[key].value) > prec * abs(param.value)  # type: ignore
                       for key, param in self.param_dict.items() if isinstance(param.value, float))  # type: ignore

    def __ne__(self, other: Block, prec: float = 1e-4):  # type: ignore
        return not self.__eq__(other, prec)

    def append(self, obj: Parameter):

        if not hasattr(self, 'name'):  # can happen if loeaded from pickle
            self.__init__(obj.lhablock)
        assert not obj.lhablock or obj.lhablock == self.name

        # The following line seems/is stupid but allow to pickle/unpickle this object
        # this is important for madspin (in gridpack mode)
        if not hasattr(self, 'param_dict'):
            self.param_dict = {}

        if tuple(obj.lhacode) in self.param_dict:
            if self.param_dict[tuple(obj.lhacode)].value != obj.value:
                raise InvalidParamCard('%s %s is already define to %s impossible to assign %s' %
                                       (self.name, obj.lhacode, self.param_dict[tuple(obj.lhacode)].value, obj.value))
            return
        list[Parameter].append(self, obj)
        # update the dictionary of key
        self.param_dict[tuple(obj.lhacode)] = obj

    def create_param_dict(self):
        """create a link between the lhacode and the Parameter"""
        for param in self:
            self.param_dict[tuple(param.lhacode)] = param

        return self.param_dict

    def def_scale(self, scale: float):
        """ """
        self.scale = scale

    def load_str(self, text: str):
        "set inforamtion from the line"

        if '#' in text:
            data, self.comment = text.split('#', 1)
        else:
            data, self.comment = text, ""

        data = data.lower()
        data = data.split()
        self.name = data[1]  # the first part of data is model
        if len(data) == 3:
            if data[2].startswith('q='):
                # the last part should be of the form Q=
                self.scale = float(data[2][2:])
            elif self.name == 'qnumbers':
                self.name += ' %s' % data[2]
        elif len(data) == 4 and data[2] == 'q=':
            # the last part should be of the form Q=
            self.scale = float(data[3])

        return self

    def keys(self):
        """returns the list of id define in this blocks"""

        return [p.lhacode for p in self]

    def __str__(self, precision: str = ''):
        """ return a str in the SLAH format """

        text = """###################################""" + \
               """\n## INFORMATION FOR %s""" % self.name.upper() +\
               """\n###################################\n"""
        # special case for decay chain
        if self.name == 'decay':
            for param in self:
                pid = param.lhacode[0]
                param.set_block('decay')
                text += str(param) + '\n'
                if pid in self.decay_table:
                    text += str(self.decay_table[pid])+'\n'
            return text
        elif self.name.startswith('decay'):
            text = ''  # avoid block definition
        # general case
        elif not self.scale:
            text += 'BLOCK %s # %s\n' % (self.name.upper(), self.comment)
        else:
            text += 'BLOCK %s Q= %e # %s\n' % (self.name.upper(),
                                               self.scale, self.comment)

        text += '\n'.join([param.__str__(precision) for param in self])
        return text + '\n'


class ParamCard(dict[str, Block]):
    """ a param Card: list of Block """
    mp_prefix = 'MP__'

    header = \
        """######################################################################\n""" + \
        """## PARAM_CARD AUTOMATICALY GENERATED BY GAMMALOOP                 ####\n""" + \
        """######################################################################\n"""

    def __init__(self, input_path: str | ParamCard | None = None):
        dict[str, Block].__init__(self, {})
        self.order: list[Block] = []
        self.not_parsed_entry: list[str] = []

        if isinstance(input_path, ParamCard):
            i_p: ParamCard = input_path
            self.read(i_p.write())
            self.input_path = input_path.input_path
        else:
            self.input_path = input_path
            if input_path:
                self.read(input_path)

    def analyze_param_card(self, model: Model | None = None) -> tuple[dict[str, list[tuple[str, tuple[int, ...]]]], dict[tuple[str, tuple[int, ...]], str]]:
        """ Analyzes the comment of the parameter in the param_card and returns
        a dictionary with parameter names in values and the tuple (lhablock, id)
        in value as well as a dictionary for restricted values.
        WARNING: THIS FUNCTION RELIES ON THE FORMATTING OF THE COMMENT IN THE
        CARD TO FETCH THE PARAMETER NAME. This is mostly ok on the *_default.dat
        but typically dangerous on the user-defined card."""

        pname2block: dict[str, list[tuple[str, tuple[int, ...]]]] = {}
        restricted_value: dict[tuple[str, tuple[int, ...]], str] = {}

        for bname, block in self.items():
            for lha_id, param in block.param_dict.items():
                all_var = []
                comment = param.comment
                if model is not None:
                    model_parameter = model.get_parameter_from_lha_specification(
                        bname, lha_id)
                    if model_parameter is not None:
                        all_var = [model_parameter.name.lower()]

                if len(all_var) == 0:
                    # treat merge parameter
                    if comment.strip().startswith('set of param :'):
                        all_var = list(re.findall(
                            r'''[^-]1\*(\w*)\b''', comment))
                    # just the variable name as comment
                    elif len(comment.split()) == 1:
                        all_var = [comment.strip().lower()]
                    # either contraction or not formatted
                    else:
                        split = comment.split()
                        if len(split) > 2 and split[1] == ':':
                            # NO VAR associated
                            restricted_value[(bname, lha_id)
                                             ] = ' '.join(split[1:])
                        elif len(split) == 2:
                            if re.search(r'''\[[A-Z]\]eV\^''', split[1]):
                                all_var = [comment.strip().lower()]
                        elif len(split) >= 2 and split[1].startswith('('):
                            all_var = [split[0].strip().lower()]
                        else:
                            if not bname.startswith('qnumbers') and comment != '':
                                logger.debug("information not recognized for %s %s : %s",
                                             bname, lha_id, comment)
                            continue

                for var in all_var:
                    var = var.lower()
                    if var in pname2block:
                        pname2block[var].append((bname, lha_id))
                    else:
                        pname2block[var] = [(bname, lha_id)]

        return (pname2block, restricted_value)

    def read(self, input_path: Any):
        """ read a card and full this object with the content of the card """

        if isinstance(input_path, str):
            if '\n' in input_path:
                input = StringIO(input_path)
            else:
                input = open(input_path)
        else:
            input = input_path  # Use for banner loading and test

        cur_block: Block | None | str = None
        for line in input:
            line = line.strip()
            if not line or line[0] == '#':
                continue
            line = line.lower()
            if line.startswith('block'):
                cur_block = Block()
                cur_block.load_str(line)
                self.append(cur_block)
                continue

            if line.startswith('decay'):
                if not self.has_block('decay'):
                    cur_block = Block('decay')
                    self.append(cur_block)
                else:
                    cur_block = self['decay']  # type: ignore
                param = Parameter()
                param.set_block(cur_block.name)  # type: ignore
                param.load_str(line[6:])
                cur_block.append(param)  # type: ignore
                continue

            if line.startswith('xsection') or cur_block == 'notparsed':
                cur_block = 'notparsed'
                self.not_parsed_entry.append(line)
                continue

            if cur_block is None:
                continue

            if cur_block.name == 'decay':  # type: ignore
                # This is a decay table
                decay_id = cur_block[-1].lhacode[0]  # type: ignore
                cur_block = Block('decay_table_%s' % decay_id)  # type: ignore
                self['decay'].decay_table[id] = cur_block  # type: ignore

            if cur_block.name.startswith('decay_table'):  # type: ignore
                param = Parameter()
                param.load_decay(line)
                try:
                    cur_block.append(param)  # type: ignore
                except InvalidParamCard:
                    pass
            else:
                param = Parameter()
                param.set_block(cur_block.name)  # type: ignore
                param.load_str(line)
                cur_block.append(param)  # type: ignore

        return self

    def __setitem__(self, name: str, value: Block) -> None:

        return dict[str, Block].__setitem__(self, name.lower(), value)

    def __getitem__(self, name: str) -> Block:
        return dict[str, Block].__getitem__(self, name.lower())

    def write(self, outpath: None | str | ParamCard = None, precision: str = '') -> str:
        """schedular for writing a card"""

        # order the block in a smart way
        blocks: list[Block] = self.order_block()
        text = self.header
        text += ''.join([block.__str__(precision) for block in blocks])
        text += '\n'
        text += '\n'.join(self.not_parsed_entry)
        if not outpath:
            return text
        elif isinstance(outpath, str):
            open(outpath, 'w').write(text)
            return ''
        else:
            return outpath.write(text)  # for test purpose

    def create_diff(self, new_card: ParamCard):
        """return a text file allowing to pass from this card to the new one
           via the set command"""

        diff = ''
        for blockname, block in self.items():
            for param in block:
                lhacode = param.lhacode
                value = param.value
                new_value = new_card[blockname].get(lhacode).value
                if isinstance(value, (float, int)) and isinstance(new_value, (float, int)):
                    same = (float(abs(value+new_value)) > 0. and not float(abs(value-new_value)/abs(
                        value+new_value)) < 1.0e-6) or (float(value) == 0. and float(new_value) == 0.)
                else:
                    same = value == new_value
                if same:
                    lhacode = ' '.join([str(i) for i in lhacode])
                    diff += 'set param_card %s %s %s # orig: %s\n' % \
                        (blockname, lhacode, new_value, value)
        return diff

    def get_value(self, blockname: str, lhecode: int | tuple[int, ...], default: float | None = None) -> float | str:
        try:
            return self[blockname].get(lhecode).value
        except KeyError:
            if blockname == 'width':
                blockname = 'decay'
                return self.get_value(blockname, lhecode, default=default)
            elif default is not None:
                return default
            raise

    def convert_to_complex_mass_scheme(self):
        """ Convert this param_card to the convention used for the complex mass scheme:
        This includes, removing the Yukawa block if present and making sure the EW input
        scheme is (MZ, MW, aewm1). """

        # The yukawa block is irrelevant for the CMS models, we must remove them
        if self.has_block('yukawa'):
            # Notice that the last parameter removed will also remove the block.
            for lhacode in [param.lhacode for param in self['yukawa']]:
                self.remove_param('yukawa', lhacode)

        # Now fix the EW input scheme
        EW_input: dict[tuple[str, tuple[int, ...]], float | None] = {('sminputs', (1,)): None,
                                                                     ('sminputs', (2,)): None,
                                                                     ('mass', (23,)): None,
                                                                     ('mass', (24,)): None}
        for block, lhaid in EW_input.keys():
            try:
                EW_input[(block, lhaid)] = self[block].get(  # type: ignore
                    lhaid).value  # type: ignore
            except:
                pass

        # Now specify the missing values. We only support the following EW
        # input scheme:
        # (alpha, GF, MZ) input
        internal_param = [key for key,
                          value in EW_input.items() if value is None]
        if len(internal_param) == 0:
            # All parameters are already set, no need for modifications
            return

        if len(internal_param) != 1:
            raise InvalidParamCard(' The specified EW inputs has more than one' +
                                   ' unknown: [%s]' % (','.join([str(elem) for elem in internal_param])))

        if not internal_param[0] in [('mass', (24,)), ('sminputs', (2,)),
                                     ('sminputs', (1,))]:
            raise InvalidParamCard(' The only EW input scheme currently supported' +
                                   ' are those with either the W mass or GF left internal.')

        # Now if the Wmass is internal, then we must change the scheme
        if internal_param[0] == ('mass', (24,)):
            aewm1 = EW_input[('sminputs', (1,))]
            Gf = EW_input[('sminputs', (2,))]
            Mz = EW_input[('mass', (23,))]
            try:
                Mw = math.sqrt((Mz**2/2.0)+math.sqrt((Mz**4/4.0)-((  # type: ignore
                              (1.0/aewm1)*math.pi*Mz**2)/(Gf*math.sqrt(2.0)))))  # type: ignore
            except Exception as _exc:
                raise InvalidParamCard(
                    'The EW inputs 1/a_ew=%{aewm1:f}, Gf=%{Gf:f}, Mz=%{Mz:f} are inconsistent')
            self.remove_param('sminputs', (2,))
            self.add_param('mass', (24,), Mw, 'MW')

    def append(self, obj: Block):
        """add an object to this"""

        assert isinstance(obj, Block)
        self[obj.name] = obj
        if not obj.name.startswith('decay_table'):
            self.order.append(obj)

    def has_block(self, name: str) -> bool:
        return name in self

    def order_block(self) -> list[Block]:
        """ reorganize the block """
        return self.order

    def rename_blocks(self, name_dict: dict[str, str]):
        """ rename the blocks """

        for old_name, new_name in name_dict.items():
            self[new_name] = self.pop(old_name)
            self[new_name].name = new_name
            for param in self[new_name]:
                param.lhablock = new_name

    def remove_block(self, name: str):
        """ remove a blocks """
        assert len(self[name]) == 0
        [self.order.pop(i) for i, b in enumerate(self.order) if b.name == name]
        self.pop(name)

    def remove_param(self, block: str, lhacode: int | tuple[int, ...]) -> None:
        """ remove a parameter """
        if self.has_param(block, lhacode):  # type: ignore
            self[block].remove(lhacode)  # type: ignore
            if len(self[block]) == 0:
                self.remove_block(block)

    def has_param(self, block: str, lhacode: int):
        """check if param exists"""

        try:
            self[block].get(lhacode)
        except:
            return False
        else:
            return True

    def add_param(self, block: str, lha: int | tuple[int, ...], value: float, comment: str = ''):

        parameter = Parameter(block=block, lhacode=lha, value=value,  # type: ignore
                              comment=comment)
        try:
            new_block = self[block]
        except KeyError:
            # If the new block didn't exist yet
            new_block = Block(block)
            self.append(new_block)
        new_block.append(parameter)


class ParamCardWriter(object):

    header = \
        """######################################################################\n""" + \
        """## PARAM_CARD AUTOMATICALY GENERATED BY GAMMALOOP                   ##\n""" + \
        """######################################################################\n"""

    sm_pdg = [1, 2, 3, 4, 5, 6, 11, 12, 13, 13, 14, 15, 16, 21, 22, 23, 24, 25]
    qnumber_template = """Block QNUMBERS %(pdg)d  # %(name)s 
        1 %(charge)d  # 3 times electric charge
        2 %(spin)d  # number of spin states (2S+1)
        3 %(color)d  # colour rep (1: singlet, 3: triplet, 8: octet)
        4 %(antipart)d  # Particle/Antiparticle distinction (0=own anti)\n"""

    @staticmethod
    def write(output_path: Path, model: Model, generic: bool = True):
        """write a valid param_card.dat"""

        parameters = [p for p in model.get_external_parameters()
                      if p.lhablock is not None and p.lhacode is not None]

        dep_mass, dep_width = ParamCardWriter.define_not_dep_param(
            parameters, model.particles, generic)

        with open(output_path, 'w') as out:
            out.write(ParamCardWriter.header)
            ParamCardWriter.write_card(
                model.particles, parameters, dep_mass, dep_width, out, generic)

    @staticmethod
    def define_not_dep_param(parameters: list[ModelParameter], particles: list[Particle], generic: bool = True) -> tuple[list[tuple[Particle, ModelParameter]], list[tuple[Particle, ModelParameter]]]:
        """define self.dep_mass and self.dep_width in case that they are 
        requested in the param_card.dat"""

        dep_mass = []
        dep_width = []
        if generic:
            dep_mass = [(part, part.mass) for part in particles
                        if part.pdg_code > 0 and
                        part.mass not in parameters]
            dep_width = [(part, part.width) for part in particles
                         if part.pdg_code > 0 and
                         part.width not in parameters]

        return dep_mass, dep_width

    @staticmethod
    def order_param(obj1: ModelParameter, obj2: ModelParameter):
        """ order parameter of a given block """

        if obj1.lhacode is None or obj2.lhacode is None:
            return 0

        maxlen = min([len(obj1.lhacode), len(obj2.lhacode)])

        for i in range(maxlen):
            if obj1.lhacode[i] < obj2.lhacode[i]:
                return -1
            elif obj1.lhacode[i] == obj2.lhacode[i]:
                return 0
            else:
                return 1
        # identical up to the first finish
        if len(obj1.lhacode) > len(obj2.lhacode):
            return 1
        elif len(obj1.lhacode) == len(obj2.lhacode):
            return 0
        else:
            return -1

    @staticmethod
    def write_card(particles: list[Particle], parameters: list[ModelParameter], dep_mass: list[tuple[Particle, ModelParameter]], dep_width: list[tuple[Particle, ModelParameter]], out: TextIO, generic: bool = True):
        # list all lhablock
        all_lhablock = set(
            [param.lhablock for param in parameters if param.lhablock is not None])

        # ordonate lhablock alphabeticaly
        all_lhablock = list(all_lhablock)
        all_lhablock.sort(key=lambda x: x.lower())
        # put at the beginning SMINPUT + MASS + DECAY
        for name in ['DECAY', 'MASS', 'SMINPUTS']:
            if name in all_lhablock:
                all_lhablock.remove(name)
                all_lhablock.insert(0, name)

        for lhablock in all_lhablock:
            ParamCardWriter.write_block(lhablock, out)
            need_writing = [param for param in parameters if
                            param.lhablock == lhablock]
            need_writing.sort(key=cmp_to_key(ParamCardWriter.order_param))
            for param in need_writing:
                ParamCardWriter.write_param(param, lhablock, out)

            if generic:
                if lhablock in ['MASS', 'DECAY']:
                    ParamCardWriter.write_dep_param_block(
                        parameters, dep_mass, dep_width, lhablock, out)

        if generic:
            ParamCardWriter.write_qnumber(particles, out)

    @staticmethod
    def write_block(name: str, out: TextIO):
        """ write a comment for a block"""

        out.writelines(
            """\n###################################""" +
            """\n## INFORMATION FOR %s""" % name.upper() +
            """\n###################################\n"""
        )
        if name.upper() != 'DECAY':
            out.write("""Block %s \n""" % name)

    @staticmethod
    def write_param(param: ModelParameter, lhablock: str, out: TextIO):
        if param.lhablock is None or param.lhacode is None or param.value is None:
            return
        lhacode = ' '.join(['%3s' % key for key in param.lhacode])
        if lhablock != 'DECAY':
            text = """  %s %e # %s \n""" % (
                lhacode, complex(param.value).real, param.name)
        else:
            text = '''DECAY %s %e \n''' % (lhacode, complex(param.value).real)
        out.write(text)

    @staticmethod
    def write_dep_param_block(parameters: list[ModelParameter], dep_mass: list[tuple[Particle, ModelParameter]], dep_width: list[tuple[Particle, ModelParameter]], lhablock: str, out: TextIO):
        text = "##  Not dependent paramater.\n"
        text += "## Those values should be edited following analytical the \n"
        text += "## analytical expression. Some generator could simply ignore \n"
        text += "## those values and use the analytical expression\n"

        if lhablock == 'MASS':
            data = dep_mass
            prefix = " "
        else:
            data = dep_width
            prefix = "DECAY "
        for part, param in data:
            if param.value is None:
                continue
            text += """%s %s %f # %s : %s \n""" % (prefix, part.pdg_code,
                                                   param.value.real, part.name,
                                                   utils.expression_to_string(param.expression))
        out.write(text)

    @staticmethod
    def write_qnumber(particles: list[Particle], out: TextIO):
        """ write qnumber """

        text = """#===========================================================\n"""
        text += """# QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE)\n"""
        text += """#===========================================================\n\n"""

        for part in particles:
            if part.pdg_code in ParamCardWriter.sm_pdg or part.pdg_code < 0:
                continue
            text += ParamCardWriter.qnumber_template % {'pdg': part.pdg_code,
                                                        'name': part.name,
                                                        'charge': 3 * part.charge,
                                                        'spin': 2 * part.spin + 1,
                                                        'color': part.color,
                                                        'antipart': part.name != part.antiname and 1 or 0}

        out.write(text)
