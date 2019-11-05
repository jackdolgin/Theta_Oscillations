# -*- coding: utf-8 -*-
# PsychoPy is required for this experiment

from psychopy import locale_setup, prefs, gui, visual, core, data, event, logging, clock, prefs
prefs.general['audioLib'] = ['pyo']
from psychopy import sound
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy.random import random, randint, normal, shuffle
import random
import os  # handy system and path functions
import sys  # to get file system encoding
import math
from itertools import chain

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
expName = 'Theta_Oscillations_script.py'
expInfo = {'participant': '', 'session': '', 'fill':''}     #creates dictionary of experiment information
dlg = gui.DlgFromDict(dictionary = expInfo, title = expName)    #creates popup at beginning of experiment that asks for participant number
if dlg.OK == False: #says, if you hit escape/click cancel when that popup appears, then don't run the experiment; if this if statement didn't exist, experiment would run regardly of whether you hit escape/click cancel
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'analysis/data/%s/%s' % (expInfo['participant'], expInfo['participant'])    #creates data file name
# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(extraInfo = expInfo, dataFileName = filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename + '.log', level = logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file; unclear what this line does




##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##-----------------------------START CODE-------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##




##--------------------------SETUP THE WINDOW----------------------------------##

win = visual.Window(
    size = (1024, 768), color = [.15, .15, .15], fullscr = False,
    allowGUI = False, monitor = 'testMonitor', useFBO = True)
# store frame rate of monitor
f_rate = win.getActualFrameRate()
expInfo['frameRate'] = f_rate


##---------------------SETUP DIFFERENT TASKS----------------------------------##

main_keys = ["left", "down", "right"]
wm_arith_keys = ["s", "S", "d", "D"]

task = int(expInfo['session'])
extra_task = int(expInfo['fill'])


##----------------------------------------------------------------------------##
##--------------------------SET UP STIMULI------------------------------------##
##----------------------------------------------------------------------------##


##-----------------------------SHAPE SIZES------------------------------------##

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

lilsize = 20
lilsize = round_up_to_even(lilsize)

def round_up_to_lilsize_multiple(f):
    return math.ceil(f * 1.0 / lilsize) * lilsize

square_size = 210
square_size = round_up_to_lilsize_multiple(square_size)

bulge_size = 1.15
cross_size = .07

if extra_task == 2:
    sides_x = 280
    sides_y = 0
    bottom_y = 0
elif extra_task == 1:
    square_mult = 528
    sides_x = int(square_mult / 2)
    sides_y = int(square_mult * math.sqrt(3) / 6)
    bottom_y = int(-square_mult * math.sqrt(3) / 4)


##-----------------------------SHAPE COLORS-----------------------------------##

shapecolor = [.05, .05, .05]
crosscolor = [1, 1, 1]
lilcolor = [-.55, -.55, -.55]

opacity = .16


##---------------------------SOUND NOISE & LOUDNESS---------------------------##

sound_clip = sound.Sound('A')
soundfiles = [os.path.join('stimuli', 'ding.wav'),
              os.path.join('stimuli', 'chord.wav')]


##-----------------------TARGET INTERVALS & DURATION SIZES--------------------##

def to_frames(t): # converts time to frames accounting for the computer's refresh rate (aka framelength); the output is in frame rates
    return int(round(t / win.monitorFramePeriod))

stagger = to_frames(.0167 * 3)
lilduration = to_frames(.0333)
bulge_duration = to_frames(.0333)
trial_start_min = to_frames(1)
trial_start_max = to_frames(1.2)
bulge_start_min = to_frames(.4)
bulge_start_max = to_frames(.8)
lilafterbulge_constant = to_frames(.3)
squareafterlil = to_frames(1)
wmarith_pause = to_frames(.5)
wmamrith_screen = 1.5
wmarith_total = wmarith_pause + to_frames(wmamrith_screen)

def while_start(k):
    global t, frameN, key
    if event.getKeys(keyList = ["escape"]):
        core.quit()
    # get current time
    t = trialClock.getTime()
    frameN += 1  # number of completed frames (so 0 is the first frame)
    key = event.getKeys(keyList = k)

def for_start():
    global frameN, continueRoutine
    frameN = -1
    continueRoutine = True

def blockdelay(): # creates a slight delay after the instruction screen and before the start of each block so the onset of the block's first trial isn't too sudden
    frameN = -1
    while frameN < to_frames(1.5):
        frameN += 1
        win.flip()

def short_on_off(r):
    s = globals()[r]
    s.setAutoDraw(True)
    while len(event.getKeys(keyList=['space'])) == 0:
        win.flip()
    s.setAutoDraw(False)

##--------------------------------CREATE TIMERS-------------------------------##

globalClock = core.MonotonicClock()  # to track the time since experiment started
trialClock = core.Clock() #unlike globalclock, gets reset each trial


##--------------------------------SET UP TRIAL COUNT--------------------------##

def round_to_multiple(f, g): #so that the num of bins can be split among even-length blocks
    return int(math.floor(f * 1.0 / g) * g)
blocks = 8
intervals = 16 # 48
intervals = round_to_multiple(intervals, blocks)
appearances = 36 #12
liltrials = intervals * appearances

percent_catch = .1      #if value gets too high script will quit out, since the step in the range function will equal 0 which isn't allowed

catchtrials = round_to_multiple((percent_catch * liltrials) / (1 - percent_catch), 6) # 6 is because of a 2x3 during catch trials- invalid locations x absent locations
total_trials = liltrials + catchtrials
trialsperblock = total_trials / blocks

validity = .7
num_squares = 4 - extra_task
num_square_multiple = round_to_multiple(validity * appearances, num_squares) # number of times we can guarantee every interval is valid at each square
extra_valids = int((validity * liltrials - (num_square_multiple * intervals)) / num_squares) # for our desired percent of valid trials, some `appearances` had to include both valid and invalid trials-and then we randomized the intervals assigned to these two 'mixed' appearances
extra_invalids = intervals - extra_valids
wmarith_freq = .3
wmarith_trials = int(total_trials * wmarith_freq)
no_wmarith_trials = total_trials - wmarith_trials

ptrials = 20
qpractice = 15 # target detection trials that are purely practice and don't contribute to the staircasing
qtrials = 30
pcutoff = 70
qcutoff = 35
startThresh = .85
acc_aim = .65
running_staircase_length = 16


##---------------------------TRIAL MATRIX & RANDOMIZATION---------------------##

d = {}
key_attr = (("x_coords", -1, 1, 0), ("y_coords", sides_y, sides_y, bottom_y), ("all_squares", "left", "right", "bottom"))
for i in range(len(key_attr)):
    one_attr = key_attr[i]
    d[one_attr[0]] = one_attr[1:len(one_attr) + 1 - extra_task]
    # d[one_attr[0]] = one_attr[1:len(one_attr)]

def most_inval(which_attr, r):
    attr = d[which_attr] # dictionary taking either an `x` or `y` input spits out the list of possible x or y coordinates (e.g. [-1, 0, 1] for `x`)
    def toss_square(o): # creates invalid cues/bulges
        direction = -1 if which_attr == "all_squares" else 1 # flips presented order of absent squares so they are in 'antiphase' with bulge at that location
        return list((attr[:o] + attr[o + 1:])[::direction])
    extra_list = []
    invalid_list = []
    hmm = ((appearances - num_square_multiple - num_squares) / (num_squares))

    if which_attr == "all_squares":
        if r == "partially_inval":
            for i in list(range(len(attr))):
                extra_inval = toss_square(i)
                extra_list.append(np.asarray(extra_inval * (extra_valids / (len(extra_inval))) + extra_inval * (extra_invalids / len(extra_inval)))) # repeat each combo of invalid x or y coordinates to fill up the valid and invalid trials in one 'mixed' interval-sets; loop three times for the three interval-sets)
                filler = ((appearances - num_square_multiple - num_squares) / num_squares) * intervals / len(extra_inval)
                for k in list(range(len(attr))):
                    invalid_list.append(toss_square(k) * (filler / len(attr)))
        elif r == "all_valid":
            for i in list(range(len(attr))):
                extra_inval = toss_square(i)
                filler = intervals * num_square_multiple / (num_squares * len(extra_inval))
                invalid_list.append(extra_inval * filler)
    elif which_attr == "x_coords" or which_attr == "y_coords":
        for i in list(range(len(attr))):
            extra_inval = toss_square(i)
            extra_list.append(np.asarray([attr[i]] * (extra_valids / (len(extra_inval))) + extra_inval * (extra_invalids / len(extra_inval))))
            ganoosh = (appearances - num_square_multiple - num_squares) / num_squares
            filler = ganoosh * intervals / len(extra_inval)
            for k in list(range(i * ganoosh, (i + 1) * ganoosh)):
                invalid_list.append(toss_square(k % num_squares) * (filler / ganoosh))
    return np.concatenate(extra_list + invalid_list)

def peat_intervals(p, q):
    return np.repeat(p, intervals * q) # repeat each location of the cue/target (`p`) `q` interval sets

def peat_catch(p, q):
    return np.repeat(p * q, int(catchtrials / (num_squares * q))) # repeat catch trials for `p` number of trials
num_square_quotient = num_square_multiple / num_squares
target_x = np.concatenate((peat_intervals(np.concatenate((np.repeat(d["x_coords"], num_square_quotient), np.tile(d["x_coords"], (appearances - num_square_multiple) / num_squares))), 1), peat_catch(d["x_coords"], 1)))
bulge_x = np.concatenate((peat_intervals(d["x_coords"], num_square_quotient), most_inval("x_coords", "partially_inval"), peat_catch(d["x_coords"], 1)))
bulge_y = np.concatenate((peat_intervals(d["y_coords"], num_square_quotient), most_inval("y_coords", "partially_inval"), peat_catch(d["y_coords"], 1)))

target_y = np.concatenate((peat_intervals(np.concatenate((np.repeat(d["y_coords"], num_square_quotient), np.tile(d["y_coords"], (appearances - num_square_multiple) / num_squares))), 1), peat_catch(d["y_coords"], 1)))

absent_list = np.concatenate((most_inval("all_squares", "all_valid"), most_inval("all_squares", "partially_inval"), peat_catch(d["all_squares"][::-1], 2)))
intervals_range = list(range(0, intervals * stagger, stagger))
lil_timing = list(chain.from_iterable([random.sample(intervals_range, len(intervals_range)) for i in list(range(appearances))])) # randomizes order of CTI's for each set of 48; so all 48 appear before the next set, but the order for each 48 is random from set to set

if catchtrials > intervals:  # spaces out when the lilsquare comes on after the bulge for catch trials
    catch_timing = [round(i * intervals * 1.0 / (catchtrials - intervals)) for i in list(range(0, (catchtrials - intervals)))]
    catch_timing += list(range(0, intervals * stagger, stagger))
else:
    catch_timing = [round(i * intervals * 1.0 / (catchtrials)) for i in list(range(0, (catchtrials)))]
np.random.shuffle(catch_timing)

def my_randomize():
    return np.asarray([random.randrange(lilsize - square_size, square_size - lilsize) / 2 for i in list(range(total_trials))]) # randomize either x or y coordinates for target within the square (gets added to target_y for that trial)


##---------------------------WM/ARITHMETIC SETUP------------------------------##

sumlist = []
corrlist = []
def createsums(s):
    numone = random.randrange(1, 9)
    numtwo = random.randrange(1, 9)
    return sumlist.append(str(numone) + " + " + str(numtwo) + " = " + str(numone + numtwo + s))

inc_arith = [-2, -1, 1, 2]
for i in inc_arith:
    for j in range(int(total_trials / (2 * len(inc_arith)))):
        createsums(i)
        corrlist.append(False)

for j in range(total_trials - len(sumlist)):
    createsums(0)
    corrlist.append(True)

def check_correct(q):
    global acc, soundp, continueRoutine
    if q:
        acc = 1
        soundp = soundfiles[0]
    else:
        acc = 0
        soundp = soundfiles[1]
    continueRoutine = False

combined = list(zip(sumlist, corrlist))
random.shuffle(combined)
sumlist[:], corrlist[:] = zip(*combined)
wmarith = wmarith_trials * [True] + no_wmarith_trials * [False]
np.random.shuffle(wmarith)


##-----------------------CREATE EXPERIMENT MATRIX-----------------------------##

expmatrix = [target_x, # side of screen of lil
            bulge_x, # side of screen of bulge
            bulge_y, # y-value of bulge
            my_randomize(), # randomize target's x value
            my_randomize(), # randomize target's y value
            [opacity] * liltrials + [0] * catchtrials, #opacity of lil
            lil_timing + catch_timing, # CTI interval for each trial
            target_y, # y value of lil,
            absent_list, # absent big square (unless task == 3)
            wmarith, # whether wm/arithmetic present
            sumlist, # arithmetic questions (if they're called)
            corrlist] # whether the answer to the arithmetic question (if it's called) is correct

#randomization sequence
randomseq = list(range(int(total_trials)))
np.random.shuffle(randomseq)
def big_opacity(s):
    return globals()["square_{0}".format(square_absent)].setOpacity(s)

def run_wm_arith():
    return (expmatrix[9][randomseq[trial]] or rep == 0) and blocktrialcount > 1 # or rep == 0 means that for any trial in rep == 0 (practice trials) there will be a wm/arith task on the screen; blocktrialcount > 1 means it can't be the first trial of the block

def save_data(r, s):
    for i in s:
        if not r:
            resp = i[1]
        elif run_wm_arith():
            resp = globals()[i[1]]
        else:
            resp = "NA"
        thisExp.addData(i[0], resp)


##----------------------MAKE THE SEARCH ARRAY ITEMS---------------------------##

cross = visual.ShapeStim(
    win = win, vertices = [[0, 0], [0, .2306 * cross_size], [0, 0],
    [0, -.2306 * cross_size], [0, 0], [.1501 * cross_size, 0],
    [0, 0], [-.1501 * cross_size, 0], [0, 0]], lineWidth = 2,
    lineColor = crosscolor, depth = -1.0)

def create_square(s, p, c):
    return visual.Rect(win = win, units = 'pix', height = s, width = s,
        pos = p, lineColor = shapecolor, fillColor = c, depth = -1.0)

square_left = create_square(square_size, (-sides_x, sides_y), shapecolor)

square_bottom = create_square(square_size, (0, bottom_y), shapecolor)

square_right = create_square(square_size, (sides_x, sides_y), shapecolor)

bulge = create_square(square_size * 1.1, (0, 0), shapecolor)

lilsquare = create_square(lilsize, (0, 0), lilcolor)


##-------------------------------INSTRUCTION SCREEN---------------------------##

if extra_task == 2:
    num_to_word = "three"
    inst4_help1 = ""
    inst4_help2 = ""
    inst68_help = main_keys[2] + ", " + main_keys[1] +", or " + main_keys[0]

elif extra_task == 1:
    num_to_word = "two"
    inst4_help1 = "in three possible locations"
    inst4_help2 = " Those two locations will vary from trial to trial."
    inst68_help = main_keys[2] + " or " + main_keys[0]

left_inst = -8

def create_inst(x, t):
    return visual.TextStim(win = win, text = t, units = 'deg', pos = (x, 0),
        height = 1, wrapWidth = 18, fontFiles = ['Lato-Reg.ttf'])

def continue_goback(s):
    return "\n\nPress space to " + s + " or \"B\" to go back."

if extra_task == 2:
    inst1 = create_inst(0, "Welcome to the study! You will be completing a sort of target-detection task, meaning quite literally you'll be looking for a target on the screen and, when you see it, detecting it by pressing a key.\n\nPress space to continue.")

    inst2 = create_inst(left_inst, "Every trial will follow the same pattern. You'll see two squares on the screen, and when there is a target that target will always be in one of those squares (as opposed to somewhere else on the screen)" + continue_goback("continue"))

    inst3 = create_inst(left_inst, "To give you a clue about which square will host the target, one of the two squares will bulge." + continue_goback("continue"))

    inst4 = create_inst(left_inst, "Then" u"—".encode('utf-8').decode('utf-8') + "if there is a target on that trial" + u"—".encode('utf-8').decode('utf-8') + "it will briefly appear in one of the two squares. It will look like a small black flicker." + continue_goback("continue"))

    inst5 = create_inst(left_inst, "All you'll have to do is press either the left or right arrow key to indicate which square the target was in. There are going to be trials where no target appears, and so if you don't see a target during a trial, don't press anything. The trial will then end on its own." + continue_goback("continue"))

    inst6 = create_inst(0, "The task moves pretty quickly (only one second to respond after the target appears), so to get the hang of it you'll go through the following practice trials. We also ask that you keep your eye focus on the cross in the center of the screen. We ask this since part of what we're trying to study is when there are no eye-movements present. So just use your peripheral vision to detect the targets (this is called covert attention)." + continue_goback("begin"))

    main_prev = create_inst(0, "The forthcoming trials will work just like the practice trials: two big squares, a bulge, and then one second to press the arrow for the side where the target appeared, or not to push a button in case you did not see a target. Again, please use only peripheral vision to detect the target and stare only at the cross in the middle." + continue_goback("begin"))

elif extra_task == 1:
    if task == 4:
        wm_arith_task = "memory"
        a_or_an = "a"
        inst7_text = "Memory trials will ask you whether the absent location in the last two trials (where there was no square on the screen) was the same location (press \"" + wm_arith_keys[1] + "\") or whether it differed (press \"" + wm_arith_keys[3] + "\")." + continue_goback("continue")
    elif task == 2 or task == 3:
        wm_arith_task = "arithmetic"
        a_or_an = "an"
        inst7_text = "Arithmetic trials will ask you whether the total on the left is the same as (press \"" + wm_arith_keys[1] + "\") the right or whether the totals differ (press \"" + wm_arith_keys[3] + "\")." + continue_goback("continue")

    inst1 = create_inst(0, "This study features two tasks, a visual attention task and one that tests " + wm_arith_task + ". The majority of trials will be visual attention tasks, interspersed with " + wm_arith_task + " trials some of the time.\n\nPress space to continue.")

    inst2 = create_inst(left_inst, "During all visual attention trials there will be a fixation point in the center of the screen, as well as other relevant stimuli in other parts of the screen. This is a covert attention task, so we ask that you keep your eyes fixated on the cross the whole time. You can attend to the other stimuli with your (covert) attention without directing your eyes at them and keeping them fixated on the cross." + continue_goback("continue"))

    inst3 = create_inst(0, "The attention task will have three phases." + continue_goback("continue"))

    inst4 = create_inst(left_inst, "In the first phase of trials, " + num_to_word + " squares " + inst4_help1 + " will appear." + inst4_help2 + continue_goback("continue"))

    inst5 = create_inst(left_inst, "Next, a cue that bulges around one of the squares will appear and indicate with " + str(int(validity * 100)) + "% likelihood which square the target (the next phase) will appear in." + continue_goback("continue"))

    inst6 = create_inst(left_inst, "Finally, a target will appear on most trials, in the form of a small dark square, inside one of the " + num_to_word + " squares with a " + str(int(validity * 100)) + "% chance in the cued location. When it appears, indicate which square it was in by pressing the " + inst68_help + u" arrow key for the respective square. If you do not see a target, indicate this by not pressing any button." + continue_goback("continue"))

    inst7 = create_inst(left_inst, inst7_text)

    inst8 = create_inst(0, u"To get the hang of it, you’ll start with two sets of practice trials. Just for the first set, there will also be ".encode('utf-8').decode('utf-8') + a_or_an + " " + wm_arith_task + " trial following every visual attention trial.\n\nAs a reminder, indicate which square the target was in by pressing the " + inst68_help + " arrow key for the respective square" + continue_goback("begin"))

    wm_probe_text = "Press \'" + wm_arith_keys[1] + "\' for same location or \'" + wm_arith_keys[3] + "\' for different location."

    wmarith_probe = visual.TextStim(
        win = win, text = "Filler", units = 'deg', height = 1, wrapWidth = 24)

    secondpractice = create_inst(0, "In this second set of practice trials, " + wm_arith_task + " trials will occur less regularly and at a more similar rate as the rest of the experiment. Press space to begin.")

    main_prev = create_inst(0, "The forthcoming trials will work just like the practice trials: visual attention trials mostly with a target present, and then intermittent " + wm_arith_task + " trials. Again, please use only peripheral vision to view stimuli during the visual attention task and stare only at the cross in the middle." + continue_goback("begin"))

plzcexp = create_inst(0, "Please see the experimenter.")

welcmmain = create_inst(0, "Welcome to the beginning of the main experiment. This experiment will last about 50 minutes. It will feature trials split among " + str(blocks - 1) + u" breaks (which will be self-timed, so you can break as long as you’d like).\n\nPress space to continue.".encode('utf-8').decode('utf-8'))

thanks = create_inst(0, "Thank you so much for your participation! Let the experimenter know that you're finished, and he'll set up the 1-minute, post-study demographic survey.")

##--------------------------PREPARE INSTRUCTION DIAGRAMS----------------------##

if extra_task == 2:
    diagrams_inst = list(range(2,6))
elif extra_task == 1:
    diagrams_inst = [2] + list(range(4, 8))

def load_diagrams(r):
    globals()["task_diagram" + str(r)] = visual.ImageStim(win = win, image = os.path.join('stimuli', "task_" + str(task) + "_extra_" + str(extra_task), 'task_diagram' + str(r) + '.png'),pos = (.54, 0), size = (.8, 1), texRes = 256)

[load_diagrams(i) for i in diagrams_inst]

def draw(j, k, m):
    globals()[(j + "{0}").format(str(k))].setAutoDraw(m)

def inst_loop(q):
    s = (inst for inst in list(range(1, 9 - 2 * (extra_task - 1))) if inst != q)
    for i in s:
        draw("inst", i, False)
        if i in diagrams_inst:
            draw("task_diagram", i, False)
    draw("inst", q, True)
    if q in diagrams_inst:
        draw("task_diagram", q, True)





##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##--------------------------START RUNNING TRIALS------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##



##----------------------------------------------------------------------------##
##--------------------------START PRACTICE TRIALS-----------------------------##
##----------------------------------------------------------------------------##


##---------------------START PRACTICE INSTRUCTIONS----------------------------##

advance = 0 # a variable that advances the instruction screen, as well as lets them go back to see a previous instruction screen

while advance < 8 - 2 * (extra_task - 1):
    if event.getKeys(keyList = ["space"]):
        advance += 1
    elif event.getKeys(keyList = ["b"]):
        if advance > 0:
            advance -= 1
    if advance == 0:
        inst_loop(1)
    elif advance == 1:
        inst_loop(2)
    elif advance == 2:
        inst_loop(3)
    elif advance == 3:
        inst_loop(4)
    elif advance == 4:
        inst_loop(5)
    elif advance == 5:
        inst_loop(6)
    elif advance == 6:
        if extra_task == 2:
            inst6.setAutoDraw(False)
        elif extra_task == 1:
            inst_loop(7)
    elif advance == 7:
        inst_loop(8)
    elif advance == 8:
        inst8.setAutoDraw(False)

    win.flip()



##----------------------------------------------------------------------------##
##--------------------------START MAIN EXPERIMENT TRIALS----------------------##
##----------------------------------------------------------------------------##


block_count = 1
q_opacity = 0
noncatch_count = 0
repstaircase = []

# the 4 (potential) reps stand for practice arithmetic rep, practice staircase rep, staircase rep, and then main experiment rep
for rep in list(range(extra_task - 1, 4)):

    rep_acc = "NA" # just a placeholder really/ arbitrary value so we can enter the while loop
    loopscompleted = 0 # no practice trial restarts at beginning of each rep

    if rep == 0:
        abbr = "p"
    elif rep == 2:
        abbr = "q"

    # the idea of this while loop is to keep repeating qwest staircase until participant is performing well enough
    while rep_acc == "NA" or ((rep == 0 or rep == 2) and rep_acc < globals()["{0}cutoff".format(abbr)]):

        if rep == 0 or rep == 2:
            if rep == 2 and loopscompleted == 0 and extra_task == 0:
                short_on_off('secondpractice') # present instructions for second practice block (aka staircasing)
            elif loopscompleted > 0:
                if rep == 0:
                    task_text = wm_arith_task
                    task_reminder = u"— ".encode('utf-8').decode('utf-8') + inst7_text
                elif rep == 2:
                    task_text = "visual attention"
                    task_reminder = ", when the target appears, indicate which square it was in by pressing the " + inst68_help + " arrow key for the respective square. If you do not see a target, indicate this by not pressing any button."
                    startThresh += .2 # increases opacity of starting opacity each they miss accuracy threshold


                ##----TELL PTCPT TO SEE EXPERIMENTER BEFORE 2ND RESTART-------##

                if loopscompleted == 2:
                    short_on_off('plzcexp')
                    loopscompleted = 1


                ##-SHOW INSTRUCTIONS AGAIN IF PTCPT HASN'T REACHED THRESHOLD--##

                tryagain = visual.TextStim(
                    win = win, text = str(rep_acc) + "% of your " + task_text + " responses were correct, which is less than the " + str(globals()["{0}cutoff".format(abbr)]) + "% threshold. As a reminder" + task_reminder, units = 'deg', height = 1, wrapWidth = 20)

                short_on_off('tryagain')

            acclist = []

        if rep == 1 or rep == 3:
            rep_acc = "break out of while after one run-through"
            if rep == 3:
                block_count = blocks
                np.random.shuffle(randomseq) # reshuffle the order of trials so that practice/staircase trials are not in the same order as experimental trials

                continueRoutineInst = True
                advance = 0 # a variable that advances the instruction screen, as well as lets them go back to see a previous instruction screen
                while continueRoutineInst:
                    if event.getKeys(keyList = ["space"]):
                        advance += 1
                    elif event.getKeys(keyList = ["b"]):
                        if advance > 0:
                            advance -= 1
                    if advance == 0:
                        main_prev.setAutoDraw(False)
                        welcmmain.setAutoDraw(True)
                    elif advance == 1:
                        welcmmain.setAutoDraw(False)
                        main_prev.setAutoDraw(True)
                    else:
                        main_prev.setAutoDraw(False)
                        continueRoutineInst = False

                    win.flip()

        if rep != 2 or loopscompleted != 0:
            blockdelay()


        for block in list(range(block_count)):
        # until the 'if block > 0' line, sets the number of trials for each run through the upcoming for loop
            if rep == 0:
                trials = list(range(ptrials))
            elif rep == 1:
                trials = list(range(qpractice))
            elif rep == 2:
                trials = data.QuestHandler(startVal = startThresh, startValSd = .23,
                    pThreshold = acc_aim, gamma = 0.05,
                    nTrials = qtrials, minVal = .01, maxVal = 4)
            elif rep == 3:
                startingtrial = block * trialsperblock
                trials = list(range(startingtrial, startingtrial + trialsperblock)) # shift trials from qwest/practice to total trials / blocks

                # set up screen between blocks showing how many blocks are left
                if block > 0:
                    break_message = visual.TextStim(
                        win = win, text = "You've reached break " + str(block) + " of " + str(block_count - 1) + ". This break is self-timed, so whenever you're ready press spacebar to continue the study.\n\nAs a reminder, 75% of little squares will be on the same side as the frame.",
                        units = 'deg', height = 1, wrapWidth = 20)

                    short_on_off('break_message')

                    blockdelay()

            blocktrialcount = 0 # count number of trials in the block to not present wm/arith probe on the first trial



            ##----------------------------------------------------------------##
            ##------------------------TRIAL FOR LOOP BEGINS-------------------##
            ##----------------------------------------------------------------##


            for trial in trials:
                for_start()
                cross.setAutoDraw(True) # cross begins on screen; is located outside of while loop since it is on screen the entire trial
                lenkeylist = 0
                overalltime = globalClock.getTime()
                blocktrialcount += 1


                ##----------------------SET TARGET & OPACITY--------------------##

                if rep == 0 or rep == 1:
                    trialopacity = opacity * startThresh
                elif rep == 2:
                    trialopacity = opacity * trial # trial here refers to the updated value in the staircase, not the actual trial number
                    trial = ptrials + qpractice
                    ptrials += 1
                elif rep == 3:
                    extracted_opacity = expmatrix[5][randomseq[trial]]
                    if extracted_opacity != 0:
                        if noncatch_count > 0 and noncatch_count % running_staircase_length == 0: # checks whether it's been 'running_staircase_length' number of experimental trials since the last resetting of opacity/staircase and therefore time to reset it
                            repstair_avg = sum(repstaircase) * 1.0 / len(repstaircase) # repstaircase is the average accuracy over the last 'running_staircase_length' experimental trials
                            q_opacity += .6 * (acc_aim - repstair_avg) # either increases or decreases the staircase by .6 opacity * the difference between the hoped-for accuracy and the average accuracy on the  last 'running_staircase_length' experimental trials
                            repstaircase = []
                        noncatch_count += 1 # counts number of experimental (non-catch) trials before re-implementing staircase/ updated opacity
                    trialopacity = extracted_opacity * q_opacity # whether no opacity * threshold

                lilsquare.setOpacity(trialopacity)


                ##-----------------SET BIG SQUARES' OPACITY-------------------##
                if task == 3:
                    square_absent = "NA"
                elif task == 2 or task == 4:
                    square_absent = expmatrix[8][randomseq[trial]]
                    if  extra_task != 2:
                        big_opacity(0)

                ##----------------SET TARGET & BULGE LOCATIONS----------------##

                lil_side = expmatrix[0][randomseq[trial]] # right now is just -1s, 0s, and 1s
                lil_height = expmatrix[7][randomseq[trial]]
                lilsquare.setPos([lil_side * sides_x + expmatrix[3][randomseq[trial]], lil_height + expmatrix[4][randomseq[trial]]])
                bulge_side = expmatrix[1][randomseq[trial]]
                bulge.setPos([bulge_side * sides_x, expmatrix[2][randomseq[trial]]])


                ##--------------SET START & DURATION OF STIMULI---------------##

                square_start = random.randrange(trial_start_min, trial_start_max + 1)
                bulge_start = square_start + random.randrange(bulge_start_min, bulge_start_max + 1)
                bulge_end = bulge_start + bulge_duration

                lilafterbulge = expmatrix[6][randomseq[trial]]
                lilstart = bulge_end + lilafterbulge_constant + lilafterbulge
                lilend = lilstart + lilduration

                squareend = lilend + squareafterlil


                ##-------------------RESET TRIAL CLOCK------------------------##

                trialClock.reset()  # clock



                ##------------------------------------------------------------##
                ##------------------WHILE LOOP BEGINS-------------------------##
                ##------------------------------------------------------------##


                while continueRoutine and frameN <= squareend:
                    while_start(main_keys)

                    ##----------------SHAPES UPDATE---------------------------##

                    if frameN == square_start:
                        square_right.tStart = t
                        square_right.frameStart = frameN
                        for a_stim in [square_right, square_left, square_bottom][0:len(d["x_coords"])]:
                            a_stim.setAutoDraw(True)
                    elif frameN == bulge_start:
                        bulge.tStart = t
                        bulge.frameStart = frameN
                        bulge.setAutoDraw(True)
                    elif frameN == bulge_end:
                        bulge.setAutoDraw(False)
                        bulge.tEnd = t
                        bulge.frameEnd = frameN
                    elif frameN == lilstart:
                        lilsquare.tStart = t
                        lilsquare.frameStart = frameN
                        lilsquare.setAutoDraw(True)
                    elif frameN == lilend:
                        lilsquare.setAutoDraw(False)
                        lilsquare.tEnd = t
                        lilsquare.frameEnd = frameN
                    elif frameN == squareend:
                        key = ['nope']
                    if len(key) > 0:
                        if frameN < lilstart:
                            lenkeylist += 1
                        else:
                            for a_stim in [square_right, square_left, square_bottom][0:len(d["x_coords"])]:
                                a_stim.setAutoDraw(False)


                            ##-------------CHECK FOR RESPONSE-----------------##

                            check_correct((key == ['nope'] and trialopacity == 0) or (key[0] == main_keys[0] and lil_side == -1) or (key[0] == main_keys[1] and lil_side == 0) or (key[0] == main_keys[2] and lil_side == 1))
                            if rep == 2:
                                trials.addResponse(acc)
                                acclist.append(acc) # creates list of qwest/staircase accuracies to determine whether participant met the cutoff for moving onto the experimental trials


                    event.clearEvents()     # Clear the previously pressed keys; too-early key presses will automatically register as incorrect


                    ##--------CHECK ALL IF COMPONENTS HAVE FINISHED-----------##

                    win.flip()
                    if not continueRoutine:  # a component has requested a forced-end of Routine
                        break



                ##----------------RESET BIG SQUARES' OPACITY------------------##

                if task != 3 and extra_task != 2:
                    big_opacity(1)


                ##----------------FINISH WHILE LOOP LEFTOVERS-----------------##

                cross.setAutoDraw(False)

                sound_clip.setSound(soundp)
                sound_clip.play() #correct and incorrect sounds are both 250 ms


                ##-------------------RECORD DATA------------------------------##

                if rep < 3:
                    thisExp.addData('Trial', -trial)
                elif rep == 3:
                    thisExp.addData('Trial', trial + 1)
                    if trialopacity != 0:
                        repstaircase.append(acc)

                save_data(False, [
                    ['ButtonPressTimeinOverallExp', overalltime],
                    ['squareStartTime', square_right.tStart],
                    ['squareStartFrame', square_right.frameStart],
                    ['flash_circleStartTime', bulge.tStart],
                    ['flash_circleStartFrame', bulge.frameStart],
                    ['flash_circleEndTime', bulge.tEnd],
                    ['flash_circleEndFrame', bulge.frameEnd],
                    ['lilsquareStartTime', lilsquare.tStart],
                    ['lilsquareStartFrame', lilsquare.frameStart],
                    ['lilsquareEndTime', lilsquare.tEnd],
                    ['lilsquareEndFrame', lilsquare.frameEnd],
                    ['lilafterflash', lilafterbulge],
                    ['ButtonPressTime', t],
                    ['Key', key[0]],
                    ['false_presses', lenkeylist],
                    ['Acc', acc],
                    ['Opacity', trialopacity],
                    ['FlashSide', bulge_side],
                    ['CorrSide', lil_side],
                    ['SquareAbsent', square_absent]
                ])



                if extra_task == 1:
                    ##------------------------------------------------------------##
                    ##--------------SET WM/ARITHMETIC TRIALS & TIMING-------------##
                    ##------------------------------------------------------------##


                    if task == 4:
                        penultimate_absent = expmatrix[8][randomseq[trial - 1]]
                        corr_wmarith_choice = (penultimate_absent == square_absent)
                    elif task == 2 or task == 3:
                        penultimate_absent = "NA"
                        corr_wmarith_choice = expmatrix[11][randomseq[trial]]

                    if run_wm_arith():
                        for_start()
                        lenkey_wmarith = 0

                        if task == 4:
                            my_inst = wm_probe_text
                            def correct_wmarith():
                                return (key[0] in wm_arith_keys[:2] and corr_wmarith_choice) or (key[0] in wm_arith_keys[2:] and corr_wmarith_choice == False)
                        elif task == 2 or task == 3:
                            my_inst = expmatrix[10][randomseq[trial]]
                            def correct_wmarith():
                                return (key[0] in wm_arith_keys[:2] and corr_wmarith_choice) or (key[0] in wm_arith_keys[2:] and corr_wmarith_choice == False)
                        wmarith_probe.setText(my_inst)


                        ##-----------------RESET TRIAL CLOCK----------------------##

                        trialClock.reset()  # clock



                        ##--------------------------------------------------------##
                        ##----------------WHILE LOOP BEGINS-----------------------##
                        ##--------------------------------------------------------##


                        while continueRoutine and frameN <= wmarith_total:
                            while_start(wm_arith_keys)

                            if frameN == wmarith_pause:
                                wmarith_probe.setAutoDraw(True)
                                wmarith_probetStart = t
                                wmarith_probeframeStart = frameN
                            elif frameN == wmarith_total:
                                 key = ['nope']
                            if len(key) > 0:
                                if frameN < wmarith_pause:
                                    lenkey_wmarith += 1

                                else:
                                    check_correct(correct_wmarith())
                                    if rep == 0 and trial > 10:
                                        acclist.append(acc) # creates list of practice accuracies which determine whether participant met the cutoff for moving onto the experimental trials


                            ##-------CHECK ALL IF COMPONENTS HAVE FINISHED--------##

                            win.flip()
                            if not continueRoutine:  # a component has requested a forced-end of Routine
                                break



                        ##--------------FINISH WHILE LOOP LEFTOVERS---------------##

                        wmarith_probe.setAutoDraw(False)

                        sound_clip.setSound(soundp)
                        sound_clip.play() #correct and incorrect sounds are both 250 ms


                    ##-------------------RECORD DATA------------------------------##

                    key_press = key[0]
                    save_data(True, [
                        ['wmarith_probeStartTime', 'wmarith_probetStart'],
                        ['wmarith_probeStartFrame', 'wmarith_probeframeStart'],
                        ['ButtonPressTime_wmarith', 't'],
                        ['Key_wmarith', 'key_press'],
                        ['false_presses_wmarith', 'lenkey_wmarith'],
                        ['Acc_wmarith', 'acc'],
                        ['PenultimateAbsent', 'penultimate_absent'],
                        ['Corr_wmarith_choice', 'corr_wmarith_choice']
                    ])

                thisExp.nextEntry()

            if rep == 0 or rep == 2:
                rep_acc = int(np.mean(acclist) * 100)
                loopscompleted += 1
                if rep == 2:
                    q_opacity = trials.mean()





short_on_off('thanks')

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename + '.csv')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
