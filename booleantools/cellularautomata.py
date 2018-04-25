from enum import Enum as _Enum
import copy as _copy

class CellularAutomata:
    """
    A class representing a basic cellular automata using a boolean function as a rule.
    """

    def __init__(self, func, initialState):
        self.func = func
        self.state = []
        for cell in initialState:
            if cell == 1:
                self.state.append(func.field.ONE)
            elif cell == 0:
                self.state.append(func.field.ZERO)
            else:
                self.state.append(cell)

        self.time = []

    def update(self):
        """
        Updates this Cellular Automata to the next state, and appends the old state to the time
        attribute.
        """
        self.time.append(_copy.deepcopy(self.state))
        out = []
        for i in range(0, len(self.state)):
            start = i - (self.func.n//2)
            end = (i + (self.func.n//2)) % len(self.state) + 1
            out.append(self.func(*_circle_slice(self.state, start, end)))

        self.state = out
        
         
    def pretty_print(self, colorMap=None):
        """
        Prints the state over time in a format that is easily visible in a Bash command line.
        """
        colorMap = colorMap if colorMap is not None else {self.func.field.ZERO: " ", self.func.field.ONE: "\033[7m \033[0m"} 
        for row in self.time:
            for cell in row:
                print(colorMap[cell], sep="", end="")

            print()

    def get_column(self, col):
        """
        Gets a column over time of the CA.

        Args:
            col (int): An integer representing the column to get

        Returns:
            list: The column over time, with index 0 being the initial state.
        """
        return [row[col] for row in self.time]

class TargetTransitions(_Enum):
        """
            An Enum representing which states under which a target CA is allowed to change state.
        """
        TARGET_TO_TARGET = 0
        TARGET_TO_NON = 1
        NON_TO_TARGET = 2
        NON_TO_NON = 3
        
        def standard_revs():
            """
            Returns the standard revision values, allowing moving from Target state to Target State,
            Nontarget state to Nontarget state, and Nontarget State to Target State.

            Returns:
                TargetTransitions: The above states.
            """
            return TargetTransitions.TARGET_TO_TARGET, TargetTransitions.NON_TO_TARGET, TargetTransitions.NON_TO_NON

class TargetedCellularAutomata(CellularAutomata):
    def __init__(self, func, initialState, targetState, *args):
        super().__init__(func, initialState)
        self.targetState = targetState
        self.reversionValues = frozenset(args)

    def update(self):
        super().update()
        oldState = self.time[-1]
        for i in range(len(oldState)):
            change = _get_change_state(self.targetState[i], oldState[i], self.state[i])
            if change not in self.reversionValues:
                self.state[i] = oldState[i]

class FiniteStateMachine:
    """
    A class representing a finite state automata for a boolean function.
    """    

    def __init__(self, func):
        self.n = func.n
        self.func = func
        self.__generateGraph()

    def __generateGraph(self):
        self.adjMatr = []
        vertices = [_dec_to_bin(i, self.n - 1) for i in range(0, 2**(self.n - 1))]
        for i in range(0, (2**(self.n-1))):
            thisRow = []
            for j in vertices:
                a = _dec_to_bin(i, self.n - 1) 
                if a[1:] == j[:len(j) - 1]:
                    thisRow.append(self.func([a[0]] + j[:]))
                else:
                    thisRow.append(None)
            self.adjMatr.append(thisRow)

def _circle_slice(lst, start, end):
    """
    Slices a list in a loop. Treats the list as a ring and slices from beginning to end,
    looping back to the beginning of the list if necessary.
    """
    if 0 <= start < end < len(lst):
        return lst[start:end]
    elif start < 0:
        return lst[start+len(lst):] + lst[:end]
    elif end >= len(lst):
        return lst[start:] + lst[0:end-len(lst)]
    elif start > end:
        return lst[start:] + lst[:end]
    elif start == end:
        return lst[start:] + lst[:end]

    print("SLICE FAILURE: ", lst, "FROM:", start, "END:", end)
    return []


def _get_change_state(target, curr, nex):
    is_curr_tar = (target == curr)
    is_nex_tar = (target == nex)

    if is_curr_tar and not is_nex_tar:
        return TargetTransitions.TARGET_TO_NON
    elif is_nex_tar and not is_curr_tar:
        return TargetTransitions.NON_TO_TARGET
    elif is_nex_tar and is_curr_tar:
        return TargetTransitions.TARGET_TO_TARGET
    else:
        return TargetTransitions.NON_TO_NON
