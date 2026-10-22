"""Finite state machines for rewriting systems (ported from the reference module)."""
class State:

    def __init__(self, name, state_machine, state_type=None, rh_rule=None):
        self.name = name
        self.transitions = {}
        self.state_machine = state_machine
        self.state_type = state_type[0]
        self.rh_rule = rh_rule

    def add_transition(self, letter, state):
        self.transitions[letter] = state

class StateMachine:

    def __init__(self, name, automaton_alphabet):
        self.name = name
        self.automaton_alphabet = automaton_alphabet
        self.states = {}
        self.add_state('start', state_type='s')

    def add_state(self, state_name, state_type=None, rh_rule=None):
        new_state = State(state_name, self, state_type, rh_rule)
        self.states[state_name] = new_state

    def __repr__(self):
        return '%s' % self.name
