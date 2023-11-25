
from enum import Enum
import re
import logging




class TargetType(Enum):
    CAS9DSB = 1
    CAS12ADSB = 2
    CAS9ABE = 3
    CAS9CBE = 4

    def length(self):
        if self is TargetType.CAS9DSB or self is TargetType.CAS9ABE or self is TargetType.CAS9CBE:
            return 23
        elif self is TargetType.CAS12ADSB:
            return 24
        else:
            raise NameError("Unknown type " + self.name)

    def editing_window(self, is_forward):
        if self is TargetType.CAS9DSB:
            if is_forward:
                return [14,19]
            else:
                return [3,9]
        if self is TargetType.CAS9ABE or self is TargetType.CAS9CBE:

            if is_forward:
                return [2, 19]
            else:
                return [3, 21]
        elif self is TargetType.CAS12ADSB:
            if is_forward:
                return [14, 23]
            else:
                return [1, 10]
        else:
            raise NameError("Unknown type " + self.name)

    def validate_sequence(self,sequence):
        if self.length() != len(sequence):
            raise NameError("Invalid length for " + self.name)
        if self == TargetType.CAS9DSB or self == TargetType.CAS9ABE or self == TargetType.CAS9CBE:
            return sequence[0:2].upper() == "CC" or sequence[-2:].upper() == "GG"
        elif self == TargetType.CAS12ADSB:
            return sequence[0:3].upper() == "TTT" or sequence[-3:].upper() == "AAA"
        else:
            raise NameError("Unknown type " + self.name)

class Target:
    __slots__ = ["target", "crispr_type", "rc_valid", ]

    def __init__(self, target_sequence, crispr_type, reverse_complement_valid=True):
        self.target = target_sequence
        self.crispr_type = crispr_type
        self.rc_valid = reverse_complement_valid

        if not self.crispr_type.validate_sequence(self.target):
            raise TypeError("Invalid sequence " + self.target + " for type " + str(self.crispr_type))

class TargetPosition:
    __slots__ = ["target", "position", "forward_orientation"]

    def __init__(self, target, position, forward_orientation):
        self.target = target
        self.position = position
        self.forward_orientation = forward_orientation


def reverse_comp(string):
    comped = [comp(c) for c in string]
    comped.reverse()
    return "".join(comped)

def comp(base):
    if base == 'A':
        return 'T'
    elif base == 'a':
        return 't'
    if base == 'C':
        return 'G'
    elif base == 'c':
        return 'g'
    if base == 'G':
        return 'C'
    elif base == 'g':
        return 'c'
    if base == 'T':
        return 'A'
    elif base == 't':
        return 'a'
    elif base > 'a':
        return 'n'
    return 'N'


class EventCigar(Enum):
    D = 0
    I = 1
    S = 2
    NONE = 3
    WT = 4
    UNKNOWN = 5

    def __str__(self):
        return str(self.name)

    @staticmethod
    def from_str(label):
        if label.upper() == 'I':
            return EventCigar.I
        elif label.upper() == 'D':
            return EventCigar.D
        elif label.upper() == 'S':
            return EventCigar.S
        elif label.upper() == 'NONE':
            return EventCigar.NONE
        elif label.upper() == 'WT':
            return EventCigar.WT
        elif label.upper() == 'UNKNOWN':
            return EventCigar.UNKNOWN
        else:
            raise TypeError("Unable to parse EventCigar symbol: " + label)

class Event:
    def __init__(self, event_cigar, event_length, position, bases=""):
        self.event_cigar = event_cigar
        self.event_length = event_length
        if self.event_length is None and not (self.event_cigar in [EventCigar.UNKNOWN, EventCigar.WT, EventCigar.NONE]):
            raise TypeError("Event length must be >= 1 for type " + self.event_cigar)
        if self.event_length is not None and self.event_length < 1:
            raise TypeError("Event length must be >= 1 for type " + self.event_cigar)

        self.position = position
        if (self.position is None or self.position < 0) and not (self.event_cigar in [EventCigar.UNKNOWN, EventCigar.WT, EventCigar.NONE]):
            raise TypeError("Position must be >= 0 ")
        self.bases = bases
        if self.bases is not None and self.bases == "":
            raise TypeError("Event bases cannot be empty")
        if self.bases is not None:
            for x in self.bases:
                if not (x in EventCaller.FASTA_BASES):
                    raise TypeError("Invalid base: " + x)
            if len(self.bases) != event_length:
                raise TypeError("Event length and bases must be equal: " + str(len(self.bases)) + " and " + str(self.event_length))

    def __str__(self):
        if self.event_length is None:
            return str(self.event_cigar.value)
        else:
            ret = str(self.event_length)
            ret += str(self.event_cigar)
            ret += "+"
            ret += str(self.position)
            if self.bases is not None:
                ret += "+" + str(self.bases)
            return ret

    def __eq__(self, other):
        if not isinstance(other, Event):
            # don't attempt to compare against unrelated types
            return NotImplemented

        return self.event_cigar == other.event_cigar and \
            (self.event_length is None and other.event_length is None) or self.event_length == other.event_length and \
            (self.position is None and other.position is None) or self.position == other.position and \
            (self.bases is None and other.bases is None) or self.bases == other.bases

    @staticmethod
    def parse_event_string(event_string):
        # 10D+44_NONE_25D+76_NONE_11D+140_1I+177+T&3D+179
        if "_" in event_string:
            raise TypeError("Individual event strings should not have a separator (_), saw one in: " + event_string)
        return [Event.parse_single_event(x) for x in event_string.split("&")]

    @staticmethod
    def parse_single_event(event_string):
        tokens = event_string.split("+")
        if len(tokens) == 3:
            type_char = tokens[0][-1]
            type_length = int(tokens[0][0:len(tokens[0]) - 1])
            event_cigar = EventCigar.from_str(type_char)

            if event_cigar == EventCigar.I:
                return Event(event_cigar, type_length, int(tokens[1]), tokens[2])
            if event_cigar == EventCigar.S:
                return Event(event_cigar, type_length, int(tokens[1]), tokens[2])
            raise TypeError("unable to parse a INS or SCAR from a length 3 event string: " + event_string)
        if len(tokens) == 2:
            type_char = tokens[0][-1]
            type_length = int(tokens[0][0:len(tokens[0]) - 1])
            event_cigar = EventCigar.from_str(type_char)

            if event_cigar == EventCigar.D:
                return Event(event_cigar, type_length, int(tokens[1]), None)
            raise TypeError("unable to parse a DEL from a length 2 event string: " + event_string)
        if len(tokens) == 1:
            event_cigar = EventCigar.from_str(tokens[0])

            if event_cigar == EventCigar.UNKNOWN or event_cigar == EventCigar.WT or event_cigar == EventCigar.NONE:
                return Event(event_cigar, None, -1, None)
            raise TypeError("unable to parse a event from a length 1 event string: " + event_string)
        raise TypeError("unable to parse event string: " + event_string)

class EventCaller:
    """
    A class to call differences between a reference and read sequence.

    Attributes:
    ----------
    reference : str
        A string representing the reference DNA sequence.
    FASTA_BASES : list
        A class attribute that lists the valid bases in a FASTA format sequence.

    Methods:
    -------
    __init__(self, reference)
        Initializes the ReferenceDifference instance with a reference sequence.

    reference_number_masking(reference)
        Static method that masks all characters in a reference sequence that are not valid FASTA bases with 'N'.

    convert_change(self, position, base)
        Converts a single nucleotide change into a numerical representation based on its position and base.

    create_reference_difference_array(self, read, cigar)
        Compares the reference and read sequences, aligned by a given CIGAR string, and returns an array of differences.
    """

    FASTA_BASES = ["A", "C", "G", "T", "U", "I", "R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "N", "-"]

    def __init__(self, reference, targets):
        """
        Initializes the ReferenceDifference instance with the reference sequence provided in upper case.

        Parameters:
        ----------
        reference : str
            The reference sequence to be used for comparison.

        reference : str
            A list of target sequences for event calling
        """
        self.reference_original = reference
        self.reference = reference.upper()
        self.targets = targets

        # setup self.target_locations, a map of target to it's positions
        self.validate_and_discover_targets()

    def validate_and_discover_targets(self):
        target_positions = {}
        for target in self.targets:
            target_position = [TargetPosition(target, m.start(), True) for m in re.finditer(target.target.upper(), self.reference.upper())]
            target_positions[target] = target_position
            if target.rc_valid:
                target_position = [TargetPosition(target, m.start(), True) for m in
                                   re.finditer(target.target.upper(), reverse_comp(self.reference.upper()))]
                target_positions[target] = target_positions[target] + target_position
        self.target_locations = target_positions


    def call_events(self, clique_read):
        """
        Creates an array of numerical representations of the differences between the reference
        and read sequences, based on a given CIGAR string.

        Parameters:
        ----------
        clique_read : str
            The clique read object to call against the reference


        Returns:
        -------
        list
            A list of lists, where each entry represents the calls within a target site (which can be a list of calls
            for complex events)

        Raises:
        ------
        AssertionError
            If the lengths of the aligned reference and read sequences do not match.
        """
        reference_sequence = self.reference
        read_sequence = clique_read.read

        # Use regular expressions to find all matches of the cigar patterns, which corresponds to
        # numbers followed by operations in the CIGAR string.
        components = re.findall(r'\d+[A-Z\=]', clique_read.cigar)

        # Initialize lists to store the aligned sequences.
        aligned_reference_sequence = []
        aligned_read = []

        # store the outcome over each target
        events_per_target = [[] for x in self.targets]

        # Initialize pointers for seq1 and seq2 to keep track of our position in each sequence.
        reference_index, read_index = 0, 0

        # Iterate over each component found in the CIGAR string.
        for component in components:
            # Split the component into the length of the operation (a number) and the operation itself (M, I, or D).
            length, operation = int(component[:-1]), component[-1]

            if operation == 'M' or operation == '=' or operation == 'X':  # match mismatch groups.
                # Append the next 'length' number of bases from each sequence to the aligned sequences.
                aligned_reference_sequence.append(reference_sequence[reference_index:reference_index + length])
                aligned_read.append(read_sequence[read_index:read_index + length])
                # Increment the pointers by the length of the match.
                reference_index += length
                read_index += length

            elif operation == 'I' or operation == 'S':  # 'I' indicates an insertion in the alignment.
                # Append a gap of 'length' to seq1 and the corresponding bases from seq2.
                aligned_reference_sequence.append('-' * length)
                aligned_read.append(read_sequence[read_index:read_index + length])
                # Only increment the pointer for seq2 since seq1 has a gap.
                read_index += length

            elif operation == 'D' or operation == 'N':  # 'D' indicates a deletion in the alignment.
                # Append the corresponding bases from seq1 and a gap of 'length' to seq2.
                aligned_reference_sequence.append(reference_sequence[reference_index:reference_index + length])
                aligned_read.append('-' * length)
                # Only increment the pointer for seq1 since seq2 has a gap.
                reference_index += length
            else:
                raise NotImplemented("We dont have CIGAR implementation for: " + component)

    def overlapping_targets(self, event_start, event_stop):
        for target,target_positions in self.target_locations.items():
            for target_position in target_positions:
                window = target.crispr_type.editing_window(target_position.forward_orientation)
                target_editing_start = target_position.position + window[0]
                target_editing_stop: object = target_position.position + window[1]

                if (event_start <= target_editing_start <= event_stop) or \
                    (target_editing_start <= event_start <= target_editing_stop) or \
                    (event_start >= target_editing_start and event_stop <= target_editing_stop) or \
                    (event_start <= target_editing_start and event_stop >= target_editing_stop):
                    return True
        return False
