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
expInfo = {'participant': '', 'session':''}     #creates dictionary of experiment information
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



# Setup different tasks
all_keys = ["l","L","a","A", "b", "B"]

task = int(expInfo['session'])

# Setup the Window
win = visual.Window(
    size = (1024, 768), color = [.15, .15, .15], fullscr = False,
    allowGUI = False, monitor = 'testMonitor', useFBO = True)
# store frame rate of monitor
f_rate = win.getActualFrameRate()
expInfo['frameRate'] = f_rate



##----------------------------------------------------------------------------##
##--------------------------SET UP STIMULI------------------------------------##
##----------------------------------------------------------------------------##


##-----------------------------SHAPE SIZES------------------------------------##

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

lilsize = 20
lilsize = round_up_to_even(lilsize)

def round_up_to_lilsize_multiple(f):
    return math.ceil(f * 1.0 /lilsize) * lilsize

square_size = 210
square_size = round_up_to_lilsize_multiple(square_size)

bulge_size = 1.15
cross_size = .07

square_mult = 528
sides_x = int(square_mult / 2)
sides_y = int(square_mult * math.sqrt(3) / 6)
bottom_y = int(-square_mult * math.sqrt(3) / 4)


##-----------------------------SHAPE COLORS-----------------------------------##

shapecolor = [.05, .05, .05]
crosscolor = [1, 1, 1]
lilcolor = [-.35, -.35, -.35]

opacity = .11


##---------------------------SOUND NOISE & LOUDNESS---------------------------##

sound_clip = sound.Sound('A')
soundfiles = [os.path.join('stimuli','ding.wav'),
              os.path.join('stimuli','chord.wav')]


##-----------------------TARGET INTERVALS & DURATION SIZES--------------------##

def to_frames(t): # converts time to frames accounting for the computer's refresh rate (aka framelength); the output is in frame rates
    return int(round(t / win.monitorFramePeriod))

stagger = to_frames(.0167)
blockdelay = to_frames(1.5) # creates a slight delay after the instruction screen and before the start of each block so the onset of the block's first trial isn't too sudden
lilduration = to_frames(.0333)
bulge_duration = to_frames(.0333)
trial_start_min = to_frames(1)
trial_start_max = to_frames(1.2)
bulge_start_min = to_frames(.4)
bulge_start_max = to_frames(.8)
lilafterbulge_constant = to_frames(.3)
squareafterlil = to_frames(1)
wmarith_pause = to_frames(.5)
wmarith_total = wmarith_pause + to_frames(1.5)

def while_start():
    global t, frameN, key
    if event.getKeys(keyList = ["escape"]):
        core.quit()
    # get current time
    t = trialClock.getTime()
    frameN += 1  # number of completed frames (so 0 is the first frame)
    key = event.getKeys(keyList = all_keys)

def for_start():
    global frameN, continueRoutine
    frameN = -1
    continueRoutine = True


##--------------------------------CREATE TIMERS-------------------------------##

globalClock = core.MonotonicClock()  # to track the time since experiment started
trialClock = core.Clock() #unlike globalclock, gets reset each trial


##--------------------------------SET UP TRIAL COUNT--------------------------##

def round_to_multiple(f, g): #so that the num of bins can be split among even-length blocks
    return math.ceil(f * 1.0 / g) * g
blocks = 1  #this line is a placeholder; change blocksreal to change # of blocks that appear
blocksreal = 8
intervals = 48
intervals = round_to_multiple(intervals, blocksreal)
reps = 12
liltrials = intervals * reps

percent_catch = .1      #if value gets too high script will quit out, since the step in the range function will equal 0 which isn't allowed

catchtrials = round_to_multiple((percent_catch * liltrials) / (1 - percent_catch), 6) # 6 is because of a 2x3 during catch trials- invalid locations x absent locations
trials = liltrials + catchtrials
trialsperblock = trials / blocksreal

validity = .7
extra_valids = int(validity * liltrials - ((reps - (6)) * intervals)) / 3 # for our desired percent of valid trials, some `reps` had to include both valid and invalid trials-and then we randomized the 48 intervals assigned to these two 'mixed' reps=
extra_invalids = intervals - extra_valids

wmarith_freq = .2
wmarith_trials = int(trials * wmarith_freq)
no_wmarith_trials = trials - wmarith_trials

ptrials = 10
qtrials = 30
qcutoff = 35
startThresh = .85
acc_aim = .65


##---------------------------TRIAL MATRIX & RANDOMIZATION---------------------##

def most_inv(z, r):
    my_dict = d[z + "_list"] # dictionary taking either an `x` or `y` input spits out the list of possible x or y coordinates
    def toss_square(): # creates invalid cues/bulges
        s = -1 if z == "o" else 1 # flips presented order of absent squares so they are in 'antiphase' with bulge at that location
        return (my_dict[:i] + my_dict[i + 1:])[::s]
    extra_list = []
    invalid_list = []
    for i in list(range(len(my_dict))):
        if r == "mixed_inv":
            if z == "o": # when deciding location of absent square
                extra_fill = toss_square()
            else:
                extra_fill = [my_dict[i]]
            extra_list.append(np.asarray(extra_fill * (extra_valids / (len(extra_fill))) + toss_square() * (extra_invalids / 2))) # repeat each combo of invalid x or y coordinates to fill up the valid and invalid trials in one 'mixed' interval-sets; loop three times for the three interval-sets
            filler = intervals / 2 # same as the line above exept just filling up with invalid x or y coordinates, for the invalid interval-sets
        else: # just when filling the absent interval-sets for validly cued trials-one set per absent location per valid location (2x3)
            filler = intervals
        invalid_list.append(toss_square() * filler)
    return np.concatenate(extra_list + invalid_list)

d = {}
d["{0}_list".format("x")] = list(range(-1, 2))
d["{0}_list".format("y")] = [sides_y, bottom_y, sides_y]
d["{0}_list".format("o")] = ["left", "bottom", "right"]
d_x = d["x_list"]
d_y = d["y_list"]
d_o = d["o_list"]

def peat_intervals(x, y):
    return np.repeat(x, intervals * y) # repeat each location of the cue/target (`x`) `y` interval sets

def peat_catch(x, y):
    return np.repeat(x * y, int(catchtrials / (3 * y))) # repeat catch trials for `x` number of trials

target_x = np.concatenate([peat_intervals(list(np.repeat(d_x, 2)) + d_x + d_x, 1), peat_catch(d_x, 1)])

bulge_x = np.concatenate([peat_intervals(d_x, 2), most_inv("x", "mixed_inv"), peat_catch(d_x, 1)])

bulge_y = np.concatenate([peat_intervals(d_y, 2), most_inv("y", "mixed_inv"), peat_catch(d_y, 1)])

target_y = np.concatenate([peat_intervals(list(np.repeat(d_y, 2)) + d_y + d_y, 1), peat_catch(d_y, 1)])

absent_list = np.concatenate([most_inv("o", "val"), most_inv("o", "mixed_inv"), peat_catch(d_o[::-1], 2)])
intervals_range = list(range(0, intervals * stagger, stagger))
lil_timing = list(chain.from_iterable([random.sample(intervals_range, len(intervals_range)) for x in list(range(reps))])) # randomizes order of CTI's for each block of 48; so all 48 appear before the next rep, but the order for each 48 is random from rep to rep

if catchtrials > intervals:  # spaces out when the lilsquare comes on after the bulge for catch trials
    catch_timing = [round(x * intervals * 1.0 / (catchtrials - intervals)) for x in list(range(0, (catchtrials - intervals)))]
    catch_timing += list(range(0, intervals * stagger, stagger))
else:
    catch_timing = [round(x * intervals * 1.0 / (catchtrials)) for x in list(range(0, (catchtrials)))]
np.random.shuffle(catch_timing)

def my_randomize():
    return np.asarray([random.randrange(lilsize - square_size, square_size - lilsize) / 2 for i in list(range(trials))]) # randomize either x or y coordinates for target within the square (gets added to target_y for that trial)


##---------------------------WM/ARITHMETIC SETUP------------------------------##

sumlist = []
corrlist = []
def createsums(s):
    numone = random.randrange(1, 9)
    numtwo = random.randrange(1, 9)
    return sumlist.append(str(numone) + " + " + str(numtwo) + " = " + str(numone + numtwo + s))

inc_arith = [-2, -1, 1, 2]
for i in inc_arith:
    for j in range(int(trials / (2 * len(inc_arith)))):
        createsums(i)
        corrlist.append(False)

for j in range(trials - len(sumlist)):
    createsums(0)
    corrlist.append(True)

def check_correct(x):
    global acc, soundp
    if (x):
        acc = 1
        soundp = soundfiles[0]
    else:
        acc = 0
        soundp = soundfiles[1]

def wmarith_na(x):
    if expmatrix[9][randomseq[trial]]:
        return (globals()[x])
    else:
        return "NA"

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
randomseq = list(range(int(trials)))
np.random.shuffle(randomseq)

def big_opacity(s):
    return globals()["square_{0}".format(square_absent)].setOpacity(s)


##----------------------MAKE THE SEARCH ARRAY ITEMS---------------------------##

cross = visual.ShapeStim(
    win = win, vertices = [[0, 0], [0, .2306 * cross_size], [0, 0],
    [0, -.2306 * cross_size], [0, 0], [.1501 * cross_size, 0],
    [0, 0], [-.1501 * cross_size, 0], [0, 0]], lineWidth = 2,
    lineColor = crosscolor, depth = -1.0)

square_bottom = visual.Rect(
    win = win, units = 'pix', height = square_size, width = square_size,
    pos = (0, bottom_y), lineColor = shapecolor,
    fillColor = shapecolor, depth = -1.0)

square_right = visual.Rect(
    win = win, units = 'pix', height = square_size, width = square_size,
    pos = (sides_x, sides_y),
    lineColor = shapecolor, fillColor = shapecolor, depth = -1.0)

square_left = visual.Rect(
    win = win, units = 'pix', height = square_size, width = square_size,
    pos = (-sides_x, sides_y),
    lineColor = shapecolor, fillColor = shapecolor, depth = -1.0)

bulge = visual.Rect(
    win = win, units = 'pix', height = square_size * 1.1, width = square_size * 1.1,
    pos = (0, 0), lineColor = shapecolor,
    fillColor = shapecolor, depth = -1.0)

lilsquare = visual.Rect(
    win = win, units = 'pix', height = lilsize, width = lilsize,
    lineColor = shapecolor, fillColor = lilcolor, depth = -1.0)

task_diagram_big_squares = visual.ImageStim(
    win = win, image = os.path.join('stimuli', 'stim_presentation_big_squares_task_' + str(task) + '.png'),
    pos = (.54, 0), size = (.8, 1), texRes = 256)

task_diagram_bulge = visual.ImageStim(
    win = win, image = os.path.join('stimuli', 'stim_presentation_bulge_task_' + str(task) + '.png'),
    pos = (.54, 0), size = (.8, 1), texRes = 256)

task_diagram_lilsquare = visual.ImageStim(
    win = win, image = os.path.join('stimuli', 'stim_presentation_lilsquare_task_' + str(task) + '.png'),
    pos = (.54, 0), size = (.8, 1), texRes = 256)

task_diagram_response = visual.ImageStim(
    win = win, image = os.path.join('stimuli', 'stim_presentation_response_task_' + str(task) + '.png'),
    pos = (.54, 0), size = (.8, 1), texRes = 256)


##-------------------------------INSTRUCTION SCREEN---------------------------##

if task == 3:
    num_to_word = "three"
    absent_instructions = ""
else:
    num_to_word = "two"
    absent_instructions =  " There are three possible locations for these big squares, and so one of the big squares will be \'absent\' per trial."
if task == 4:
    wm_arith_task = "memory"
else:
    wm_arith_task = "arithmetic"

inst0 = visual.TextStim(
    win = win, text = "This study features two tasks, a visual attention task and one that tests " + wm_arith_task + ". These " + wm_arith_task + "trials will follow some of the visual attention trials, but there will always be at least one visual attention trial between " + wm_arith_task + " trials.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst01 = visual.TextStim(
    win = win, text = "During all visual attention trials there will be a fixation point in the center of the screen, as well as other relevant stimuli in other parts of the screen. These visual trials are covert attention tasks, so we ask that you keep your eyes fixated on the cross the whole time. You can attend to the other stimuli with your (covert) attention without directing your eyes at them and keeping them fixated on the cross.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst02 = visual.TextStim(
    win = win, text = "This study will begin with two sets of practice trials. On each practice trial, like each trial in the main experiment, you will see " + num_to_word + " big squares appear on the screen.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst04 = visual.TextStim(
    win = win, text = "For example, one trial may begin with:",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst05 = visual.TextStim(
    win = win, text = "and the next with :",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst06 = visual.TextStim(
    win = win, text = "Once these " + num_to_word + " big squares appear, one of them will bulge shortly thereafter. This bulge is a cue for detecting a target.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst1 = visual.TextStim(
    win = win, text = "This study will begin with " + str(ptrials + qtrials) + " practice trials. On each practice trial, like each trial in the main experiment, you will see " + num_to_word + " big squares appear on the screen.\n\nPress space to continue.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst2 = visual.TextStim(
    win = win, text = "One of these squares will then bulge shortly thereafter.\n\nPress space to continue or \"B\" to go back.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst3 = visual.TextStim(
    win = win, text = "Then on some trials you will see a little square inside one of the " + num_to_word + " squares, and on other trials you will not see any little square. " + str(int(validity * 100)) + "% of little squares will occur at the same square as the bulge, and you are encouraged to use this information to aid your performance.\n\nPress space to continue or \"B\" to go back.",
    units='deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst4 = visual.TextStim(
    win = win, text = "When the little square appears, you will have one second to indicate whether it was on the left, bottom, or right of the screen by pressing, respectively, the \"A\", \"B\", or \"L\" keys. If you do not see a little square, indicate this by not pressing any button. Please gaze at the cross in the center of the screen the whole time and detect the little square with your peripheral vision.\n\nPress space to continue or \"B\" to go back.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst5 = visual.TextStim(
    win = win, text = "Again, these are practice trials, and if you do well enough on the practice trials, you will move on directly to the main experiment; otherwise you'll get a second chance to improve on the practice. As a reminder, when the little square appears, you will have one second to indicate whether it was on the left, bottom, or right of the screen by pressing, respectively, the \"A\", \"B\", or \"L\" keys.\n\nPress space to begin or \"B\" to go back.",
    units = 'deg', height = 1, wrapWidth = 20)

inst6 = visual.TextStim(
    win = win, text = "Please see the experimenter.",
    units = 'deg', height = 1, wrapWidth = 20)

inst7 = visual.TextStim(
    win = win, text = "Welcome to the beginning of the main experiment. This experiment will last about 40 minutes. It will feature trials split among " + str(blocksreal - 1) + " breaks (the breaks will be self-timed, so you can take as long as you'd like during them before proceeding to the subsequent trials).\n\nPress space to continue.",
    units = 'deg', height = 1, wrapWidth = 20)

inst8 = visual.TextStim(
    win = win, text = "The forthcoming trials will work just like the practice trials: " + num_to_word + " big squares, a bulge, and then one second to push a button to answer where the little square appeared, or not to push a button in case you did not see a little square. Again, please gaze at the cross in the center of the screen the whole time and only use peripheral vision to detect the little square.\n\nPress space to begin or \"B\" to go back.",
    units = 'deg', height = 1, wrapWidth = 20)

inst9 = visual.TextStim(
    win = win, text = "Thank you so much for your participation! Let the experimenter know that you're finished, and he'll set up the 1-minute, post-study demographic survey.",
    units = 'deg', height = 1, wrapWidth = 20)

wmarith_probe = visual.TextStim(
    win = win, text = "Filler question",
    units = 'deg', height = 1, wrapWidth = 20)

inst12 = 'Was the omitted square in the same location during the last two trials?\nPress \"T\" for True or \"F\" for False.'

# inst01 = visual.TextStim(
#     win = win, text = "This study features two tasks, a visual attention task and one that (tests recall / arithmetic). These (recall/ arithmetic) trials will follow some of the visual attention trials, but there will always be at least one visual attention trial between recall/arithmetic trials.",
#     units = 'deg', height = 1, wrapWidth = 20)
#
# inst02 = visual.TextStim(
#     win = win, text = "During all visual attention trials there will be a fixation point in the center of the screen, as well as other relevant stimuli in other parts of the screen. These visual trials are covert attention tasks, so we ask that you keep your eyes fixated on the cross the whole time. You can attend to the other stimuli with your (covert) attention without directing your eyes at them and keeping them fixated on the cross.",
#     units = 'deg', height = 1, wrapWidth = 20)
#
# #inst1 #There are three possible locations for these big squares, and so one of the big squares will be `absent` per trial.
#
# # For example, one trial may begin with : and the next with :__
#
# # Once these two big squares appear, one of them will bulge shortly thereafter. This bulge is a cue for detecting a target.
#
# #






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

while advance < 5:
    if event.getKeys(keyList = ["space"]):
        advance += 1
    elif event.getKeys(keyList = ["b"]):
        if advance > 0:
            advance -= 1
    if advance == 0:
        inst2.setAutoDraw(False)
        inst3.setAutoDraw(False)
        inst4.setAutoDraw(False)
        inst5.setAutoDraw(False)
        task_diagram_bulge.setAutoDraw(False)
        task_diagram_lilsquare.setAutoDraw(False)
        task_diagram_response.setAutoDraw(False)
        inst1.setAutoDraw(True)
        task_diagram_big_squares.setAutoDraw(True)
    elif advance == 1:
        inst1.setAutoDraw(False)
        inst3.setAutoDraw(False)
        inst4.setAutoDraw(False)
        inst5.setAutoDraw(False)
        task_diagram_big_squares.setAutoDraw(False)
        task_diagram_lilsquare.setAutoDraw(False)
        task_diagram_response.setAutoDraw(False)
        inst2.setAutoDraw(True)
        task_diagram_bulge.setAutoDraw(True)
    elif advance == 2:
        inst2.setAutoDraw(False)
        inst4.setAutoDraw(False)
        inst5.setAutoDraw(False)
        task_diagram_bulge.setAutoDraw(False)
        task_diagram_response.setAutoDraw(False)
        inst3.setAutoDraw(True)
        task_diagram_lilsquare.setAutoDraw(True)
    elif advance == 3:
        inst3.setAutoDraw(False)
        inst5.setAutoDraw(False)
        task_diagram_lilsquare.setAutoDraw(False)
        inst4.setAutoDraw(True)
        task_diagram_response.setAutoDraw(True)
    elif advance == 4:
        inst4.setAutoDraw(False)
        task_diagram_response.setAutoDraw(False)
        inst5.setAutoDraw(True)
    else:
        inst5.setAutoDraw(False)


    win.flip()


##----------------------------------------------------------------------------##
##--------------------------START MAIN EXPERIMENT TRIALS----------------------##
##----------------------------------------------------------------------------##



trials = list(range(ptrials))
q_opacity = 0
qloop = 0
noncatch_count = 0
repstaircase = []

for rep in list(range(3)):
    q_acc = 0
    if rep == 2:
        blocks = blocksreal
        np.random.shuffle(randomseq) # reshuffle the order of trials so that practice/staircase trials are not in the same order as experimental trials
    while q_acc < qcutoff:
        acclist = []
        if rep == 1:


            ##------TELL PTCPT TO SEE EXPERIMENTER BEFORE 2ND RESTART---------##

            while qloop == 2:
                inst6.setAutoDraw(True)
                if(event.getKeys(keyList = ["space"])):
                    inst6.setAutoDraw(False)
                    qloop = 1
                win.flip()


            ##---SHOW INSTRUCTIONS AGAIN IF PTCPT HASN'T REACHED THRESHOLD----##

            if qloop == 1:
                inst4 = visual.TextStim(
                    win = win, text = str(q_acc) + "% of your responses were correct, which is less than the " + str(qcutoff) + "% threshold. As a reminder, when the little square appears, you will have one second to indicate whether it was on the left or right of the screen by pressing, respectively, the \"A\" or \"L\" keys. If you do not see a little square, indicate this by not pressing any button. Finally, please detect the little square with your peripheral vision and gaze at the cross in the center of the screen the whole time.\n\nPlease press the space bar to try again.",
                    units = 'deg', height = 1, wrapWidth = 20)
                continueRoutine_q = True
                inst4.setAutoDraw(True)
                while len(event.getKeys(keyList=['space'])) == 0:
                    win.flip()
                inst4.setAutoDraw(False)

                startThresh += .2 # makes starting point for staircase easier (increases opacity) each time they score less than the practice threshold

                frameN = -1
                while frameN < blockdelay:
                    frameN += 1
                    win.flip()

            trials = data.QuestHandler(startVal = startThresh, startValSd = .23,
                pThreshold = acc_aim, gamma = 0.05,
                nTrials = qtrials, minVal = .01, maxVal = 4)
        else:
            q_acc = qcutoff + 1
            frameN = -1 # number of completed frames (so 0 is the first frame)
            while frameN < blockdelay: # creates a 1.5 second blank screen that appears after the instructions but before the first practice trial begins
                frameN += 1
                win.flip()

        for block in list(range(blocks)):
            if rep == 2:
                continueRoutineInst = True

                inst7 = visual.TextStim(
                    win = win, text = "You've reached break " + str(block) + " of " + str(blocks-1) + ". This break is self-timed, so whenever you're ready press spacebar to continue the study.\n\nAs a reminder, 75% of little squares will be on the same side as the frame.",
                    units = 'deg', height = 1, wrapWidth = 20)

                advance = 0 # a variable that advances the instruction screen, as well as lets them go back to see a previous instruction screen

                while continueRoutineInst:
                    if block == 0:
                        if(event.getKeys(keyList = ["space"])):
                            advance += 1
                        elif(event.getKeys(keyList = ["b"])):
                            if advance > 0:
                                advance -= 1
                        if (advance == 0):
                            inst8.setAutoDraw(False)
                            inst7.setAutoDraw(True)
                        elif (advance == 1):
                            inst7.setAutoDraw(False)
                            inst8.setAutoDraw(True)
                        else:
                            inst8.setAutoDraw(False)
                            continueRoutineInst = False

                    else:
                        inst7.setAutoDraw(True)
                        if (event.getKeys(keyList = ["space"])):
                            continueRoutineInst = False
                            inst7.setAutoDraw(False)

                    if continueRoutineInst:
                        win.flip()

                frameN = -1
                while frameN < blockdelay:
                    frameN += 1
                    win.flip()

                startingtrial = block * trialsperblock
                trials = list(range(startingtrial, startingtrial + trialsperblock))



            ##----------------------------------------------------------------##
            ##------------------------TRIAL FOR LOOP BEGINS-------------------##
            ##----------------------------------------------------------------##


            for trial in trials:
                for_start()
                cross.setAutoDraw(True) # cross begins on screen; is located outside of while loop since it is on screen the entire trial
                lenkeylist = 0
                overalltime = globalClock.getTime()


                ##----------------------SET TARGET & OPACITY--------------------##

                if rep == 0:
                    trialopacity = opacity * startThresh
                elif rep == 1:
                    trialopacity = opacity * trial # trial here refers to the updated value in the staircase, not the actual trial number
                    trial = ptrials
                    ptrials += 1
                else:
                    extracted_opacity = expmatrix[5][randomseq[trial]]
                    if extracted_opacity > 0:
                        if noncatch_count > 0 and noncatch_count % 16 == 0:
                            repstair_avg = sum(repstaircase) * 1.0 / len(repstaircase)
                            q_opacity += .6 * (acc_aim - repstair_avg)
                            repstaircase = []
                        noncatch_count += 1 # counts number of experimental (non-catch) trials before re-implementing staircase/ updated opacity
                    trialopacity = extracted_opacity * q_opacity # whether no opacity * threshold

                lilsquare.setOpacity(trialopacity)


                ##-----------------SET BIG SQUARES' OPACITY-------------------##
                if task == 3:
                    square_absent = "NA"
                else:
                    square_absent = expmatrix[8][randomseq[trial]]
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
                    while_start()

                    ##----------------SHAPES UPDATE---------------------------##

                    if frameN == square_start:
                        square_right.tStart = t
                        square_right.frameStart = frameN
                        square_left, square_right.setAutoDraw(True)
                        square_left.setAutoDraw(True)
                        square_bottom.setAutoDraw(True)
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
                            continueRoutine = False
                            square_right.setAutoDraw(False)
                            square_left.setAutoDraw(False)
                            square_bottom.setAutoDraw(False)


                            ##-------------CHECK FOR RESPONSE-----------------##

                            check_correct((key == ['nope'] and trialopacity == 0) or (key == ['l'] and lil_side == 1) or (key == ['a'] and lil_side == -1) or (key == ['b'] and lil_side == 0))
                            if rep == 1:
                                trials.addResponse(acc)
                            acclist.append(acc)


                    event.clearEvents()     # Clear the previously pressed keys; too-early key presses will automatically register as incorrect


                    ##--------CHECK ALL IF COMPONENTS HAVE FINISHED-----------##

                    if not continueRoutine:  # a component has requested a forced-end of Routine
                        break
                    else:
                        win.flip()


                ##----------------RESET BIG SQUARES' OPACITY------------------##

                if task != 3:
                    big_opacity(1)


                ##----------------FINISH WHILE LOOP LEFTOVERS-----------------##

                cross.setAutoDraw(False)

                sound_clip.setSound(soundp)
                sound_clip.play() #correct and incorrect sounds are both 250 ms


                ##-------------------RECORD DATA------------------------------##

                if rep < 2:
                    thisExp.addData('Trial', -trial)
                else:
                    thisExp.addData('Trial', trial + 1)
                    if trialopacity > 0:
                        repstaircase.append(acc)
                thisExp.addData('ButtonPressTimeinOverallExp', overalltime)
                thisExp.addData('squareStartTime', square_right.tStart)
                thisExp.addData('squareStartFrame', square_right.frameStart)
                thisExp.addData('flash_circleStartTime', bulge.tStart)
                thisExp.addData('flash_circleStartFrame', bulge.frameStart)
                thisExp.addData('flash_circleEndTime', bulge.tEnd)
                thisExp.addData('flash_circleEndFrame', bulge.frameEnd)
                thisExp.addData('lilsquareStartTime', lilsquare.tStart)
                thisExp.addData('lilsquareStartFrame', lilsquare.frameStart)
                thisExp.addData('lilsquareEndTime', lilsquare.tEnd)
                thisExp.addData('lilsquareEndFrame', lilsquare.frameEnd)
                thisExp.addData('lilafterflash', lilafterbulge)
                thisExp.addData('ButtonPressTime', t)
                thisExp.addData('Key', key[0])
                thisExp.addData('false_presses', lenkeylist)
                thisExp.addData('Acc', acc)
                thisExp.addData('Opacity', trialopacity)
                thisExp.addData('FlashSide', bulge_side)
                thisExp.addData('CorrSide', lil_side)
                thisExp.addData('SquareAbsent', square_absent)

                if rep == 1:
                    q_opacity = trials.mean()
                    q_acc = round(np.mean(acclist) * 100)
                    qloop += 1



                ##------------------------------------------------------------##
                ##--------------SET WM/ARITHMETIC TRIALS & TIMING-------------##
                ##------------------------------------------------------------##


                if task == 4:
                    penultimate_absent = expmatrix[8][randomseq[trial - 1]]
                    last_absent = expmatrix[8][randomseq[trial]]
                    corr_wmarith_choice = (penultimate_absent == last_absent)
                else:
                    penultimate_absent = "NA"
                    last_absent = "NA"
                    corr_wmarith_choice = (expmatrix[11][randomseq[trial]])

                if expmatrix[9][randomseq[trial]]:
                    for_start()
                    lenkey_wmarith = 0

                    if task == 4:
                        my_inst = inst12
                        correct_wmarith = ((key == "T" and corr_wmarith_choice) or (key == "F" and corr_wmarith_choice == False))
                    else:
                        my_inst = expmatrix[10][randomseq[trial]]
                        correct_wmarith = (key == "T" and expmatrix[11][randomseq[trial]]) or (key == "F" and expmatrix[11][randomseq[trial]] == False)
                    wmarith_probe.setText(my_inst)


                    ##-----------------RESET TRIAL CLOCK----------------------##

                    trialClock.reset()  # clock



                    ##--------------------------------------------------------##
                    ##----------------WHILE LOOP BEGINS-----------------------##
                    ##--------------------------------------------------------##


                    while continueRoutine and frameN <= wmarith_total:
                        while_start()

                        if frameN == wmarith_pause:
                            wmarith_probe.setAutoDraw(True)
                            wmarith_probe.tStart = t
                            wmarith_probe.frameEnd = frameN
                        elif frameN == wmarith_total:
                             key = ['nope']
                        if len(key) > 0:
                            if frameN < wmarith_total:
                                lenkey_wmarith += 1


                            else: ##------------CHECK FOR RESPONSE------------##

                                check_correct(correct_wmarith)


                        ##-------CHECK ALL IF COMPONENTS HAVE FINISHED--------##

                        if not continueRoutine:  # a component has requested a forced-end of Routine
                            break
                        else:
                            win.flip()


                    ##--------------FINISH WHILE LOOP LEFTOVERS---------------##

                    wmarith_probe.setAutoDraw(False)

                    sound_clip.setSound(soundp)
                    sound_clip.play() #correct and incorrect sounds are both 250 ms


                ##-------------------RECORD DATA------------------------------##
                # thisExp.addData('wmarith_probeStartTime', wmarith_na('wmarith_probe.tStart'))
                # thisExp.addData('wmarith_probeStartFrame', wmarith_na('wmarith_probe.frameStart'))
                thisExp.addData('ButtonPressTime_wmarith', wmarith_na('t'))
                # thisExp.addData('Key_wmarith', wmarith_na('key[0]'))
                thisExp.addData('false_presses_wmarith', wmarith_na('lenkey_wmarith'))
                thisExp.addData('Acc_wmarith', wmarith_na('acc'))
                thisExp.addData('PenultimateAbsent', wmarith_na('penultimate_absent'))
                thisExp.addData('LastAbsent', wmarith_na('last_absent'))
                thisExp.addData('Corr_wmarith_choice', wmarith_na('corr_wmarith_choice'))
                thisExp.nextEntry()





inst9.setAutoDraw(True)
while len(event.getKeys(keyList=['space'])) == 0:
    win.flip()
inst9.setAutoDraw(False)

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename + '.csv')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
