"""This module tests function reverse_words
provided by module mod.py."""
import unittest
from clique import callers
from clique.callers import TargetType, Target, reverse_comp, EventCaller, EventCigar, Event
from pprint import pprint
import logging


class TargetTypeTest(unittest.TestCase):
    def testTargetTypeCAS9DSB(self):
        target_type = TargetType.CAS9_DSB

        self.assertEqual(target_type.length(), 23)
        self.assertEqual(target_type.validate_sequence("ACGTAACGTAACGTAACGTACGG"), True)
        self.assertEqual(target_type.validate_sequence("ACGTAACGTAACGTAACGTACAT"), False)
        self.assertEqual(target_type.validate_sequence("CCGTAACGTAACGTAACGTACAT"), True)

    def testTargetTypeCAS9BaseEditors(self):
        target_type = TargetType.CAS9_ABE

        self.assertEqual(target_type.length(), 23)
        self.assertEqual(target_type.validate_sequence("ACGTAACGTAACGTAACGTACGG"), True)
        self.assertEqual(target_type.validate_sequence("ACGTAACGTAACGTAACGTACAT"), False)
        self.assertEqual(target_type.validate_sequence("CCGTAACGTAACGTAACGTACAT"), True)

        target_type = TargetType.CAS9_CBE

        self.assertEqual(target_type.length(), 23)
        self.assertEqual(target_type.validate_sequence("ACGTAACGTAACGTAACGTACGG"), True)
        self.assertEqual(target_type.validate_sequence("ACGTAACGTAACGTAACGTACAT"), False)
        self.assertEqual(target_type.validate_sequence("CCGTAACGTAACGTAACGTACAT"), True)

    def testTargetTypeCAS12A(self):
        target_type = TargetType.CAS12A_DSB

        self.assertEqual(target_type.length(), 24)
        self.assertEqual(target_type.validate_sequence("TTTACGTAACGTAACGTAACGTAC"), True)
        self.assertEqual(target_type.validate_sequence("ACGTAACGTAACGTAACGTACAAT"), False)
        self.assertEqual(target_type.validate_sequence("TTAACGTAACGTAACGTACATAAA"), True)

    def testTargetTypeCAS9_PAL_ABE(self):
        target_type = TargetType.CAS9_PAL_ABE

        self.assertEqual(target_type.length(), 26)
        self.assertEqual(target_type.validate_sequence("CCAAAAAATTTTTAAAAATTTTTCGG"), True)
        self.assertEqual(target_type.validate_sequence("CAAAAAAATTTTTAAAAATTTTTCGG"), False)
        self.assertEqual(target_type.validate_sequence("CCAAAAAATTTTTAAAAATTTTTCGA"), False)

class EventCigarTypeTest(unittest.TestCase):
    def testEventCigar(self):

        self.assertEqual(EventCigar.from_str("D"), EventCigar.D)
        self.assertEqual(EventCigar.from_str("d"), EventCigar.D)
        self.assertEqual(EventCigar.from_str("I"), EventCigar.I)
        self.assertEqual(EventCigar.from_str("i"), EventCigar.I)
        self.assertEqual(EventCigar.from_str("S"), EventCigar.S)
        self.assertEqual(EventCigar.from_str("s"), EventCigar.S)
        self.assertEqual(EventCigar.from_str("NONE"), EventCigar.NONE)
        self.assertEqual(EventCigar.from_str("NonE"), EventCigar.NONE)
        self.assertEqual(EventCigar.from_str("UNKNOWN"), EventCigar.UNKNOWN)
        self.assertEqual(EventCigar.from_str("UnKNOWN"), EventCigar.UNKNOWN)
        self.assertEqual(EventCigar.from_str("WT"), EventCigar.WT)
        self.assertEqual(EventCigar.from_str("wt"), EventCigar.WT)
        with self.assertRaises(TypeError):
            EventCigar.from_str("e")
        with self.assertRaises(TypeError):
            EventCigar.from_str("DD")
        with self.assertRaises(TypeError):
            EventCigar.from_str("II")
        with self.assertRaises(TypeError):
            EventCigar.from_str("iI")
        with self.assertRaises(TypeError):
            EventCigar.from_str("notathing")

class EventTest(unittest.TestCase):
    def testEvent(self):
        # simple events
        self.assertEqual(Event.parse_event_string("5D+100"), [Event(EventCigar.D, 5, 100, None)])
        self.assertEqual(Event.parse_event_string("4I+10+AAAA"), [Event(EventCigar.I, 4, 10, "AAAA")])
        self.assertEqual(Event.parse_event_string("5S+120+TTTTT"), [Event(EventCigar.S, 5, 120, "TTTTT")])
        self.assertEqual(Event.parse_event_string("NONE"), [Event(EventCigar.NONE, None, None, None)])
        self.assertEqual(Event.parse_event_string("UNKNOWN"), [Event(EventCigar.UNKNOWN, None, None, None)])
        self.assertEqual(Event.parse_event_string("WT"), [Event(EventCigar.WT, None, None, None)])

        # complex event strings
        self.assertEqual(Event.parse_event_string("5D+5&WT"), [Event(EventCigar.D, 5, 5, None), Event(EventCigar.WT, None, None, None)])
        self.assertEqual(Event.parse_event_string("4I+5+AAAA&UNKNOWN"), [Event(EventCigar.I, 4, 5, "AAAA"), Event(EventCigar.WT, None, None, None)])

        # some smattering of error conditions
        with self.assertRaises(TypeError):
            Event.parse_event_string("5I+100")
        with self.assertRaises(TypeError):
            Event.parse_event_string("5I+100+A")
        with self.assertRaises(TypeError):
            Event.parse_event_string("5S+100+A")
        with self.assertRaises(TypeError):
            Event.parse_event_string("5D")
        with self.assertRaises(TypeError):
            Event.parse_event_string("5S")
        with self.assertRaises(TypeError):
            Event.parse_event_string("AS")
        with self.assertRaises(TypeError):
            Event.parse_event_string("5I+100+AZ")
        with self.assertRaises(ValueError):
            Event.parse_event_string("WT+5")

class TargetTest(unittest.TestCase):
    def testTargetCreation(self):

        # some correct target creation steps
        target = Target("ACGTAACGTAACGTAACGTACGG", TargetType.CAS9_ABE)
        target = Target("TTTACGTAACGTAACGTAACGTAC", TargetType.CAS12A_DSB)
        target = Target("ACGTAACGTAACGTAACGTACGG", TargetType.CAS9_DSB)

        # assertions
        with self.assertRaises(TypeError):
            Target("ACGTAACGTAACGTAACGTACAG", TargetType.CAS9_ABE)
            Target("ACGTAACGTAACGTAACGTACAG", TargetType.CAS12A_DSB)


class ReverseCompTest(unittest.TestCase):
    def testReverseComp(self):
        # some correct target creation steps
        self.assertEqual(reverse_comp("AAAAaaaa"), "ttttTTTT")
        self.assertEqual(reverse_comp("AAAAzzzZ"), "NnnnTTTT")
        self.assertEqual(reverse_comp("ACGTAacgt"), "acgtTACGT")

class EventCallerTest(unittest.TestCase):
    def testSimpleFindInReference(self):
        target = Target("ACGTAACGTAACGTAACGTACGG", TargetType.CAS9_DSB)
        eventCaller = EventCaller("ACGTAACGTAACGTAACGTACGGAAAACGTAACGTAACGTAACGTACGGAAA",[target])
        self.assertEqual(len(eventCaller.target_locations), 1,
                         msg='targets = {0}'.format(pprint(eventCaller.target_locations)))
        self.assertEqual(len(eventCaller.target_locations[target]), 2,
                         msg='targets = {0}'.format(pprint(eventCaller.target_locations)))


    def testMixFindInReference(self):
        target1 = Target("ACGTAACGTAACGTAACGTACGG", TargetType.CAS9_DSB)
        target2 = Target("ACGTAACGTAACGTAACGTTAAAA", TargetType.CAS12A_DSB)

        # ACGTAACGTAACGTAACGTACGGACGTAACGTAACGTAACGTTAAAAAAAACGTAACGTAACGTAACGTACGGAAACCGTACGTTACGTTACGTTACGT
        # ACGTAACGTAACGTAACGTACGG                           ACGTAACGTAACGTAACGTACGG
        #                        ACGTAACGTAACGTAACGTTAAAA
        #          |10       |20       |30       |40       |50       |60

        eventCaller = EventCaller("ACGTAACGTAACGTAACGTACGGACGTAACGTAACGTAACGTTAAAAAAAACGTAACGTAACGTAACGTACGGAAACCGTACGTTACGTTACGTTACGT",[target1,target2])
        self.assertEqual(len(eventCaller.target_locations), 2,
                         msg='targets = {0}'.format(pprint(eventCaller.target_locations)))
        self.assertEqual(len(eventCaller.target_locations[target1]), 3,
                         msg='targets = {0}'.format(pprint(eventCaller.target_locations)))
        self.assertEqual(len(eventCaller.target_locations[target2]), 1,
                         msg='targets = {0}'.format(pprint(eventCaller.target_locations)))


    def testOverlappingTargets(self):
        #      target 1                   target 2
        # [ 0 ---------------- 22][ 23 ---------------- 45]
        #
        # ACGTAACGTAACGTAACGTACGGACGTAACGTAACGTAACGTTAAAAAAAACGTAACGTAACGTAACGTACGGAAACCGTACGTTACGTTACGTTACGT
        #          |10       |20       |30       |40        |50        |60       |70       |80
        target1 = Target("ACGTAACGTAACGTAACGTACGG", TargetType.CAS9_DSB) # position 0, 50
        target2 = Target("ACGTAACGTAACGTAACGTTAAAA", TargetType.CAS12A_DSB) # 23



        eventCaller = EventCaller("ACGTAACGTAACGTAACGTACGGACGTAACGTAACGTAACGTTAAAAAAAACGTAACGTAACGTAACGTACGGAAACCGTACGTTACGTTACGTTACGT",[target1,target2])
        self.assertEqual(eventCaller.overlapping_targets(10,20), True)
        self.assertEqual(eventCaller.overlapping_targets(30,40), True)
        self.assertEqual(eventCaller.overlapping_targets(47, 49), False)
        self.assertEqual(eventCaller.overlapping_targets(45, 50), True)
        self.assertEqual(eventCaller.overlapping_targets(75, 80), False)