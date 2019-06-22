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
expName = 'Theta_Oscillations_Exp_1.py'
expInfo = {'participant': ''}     #creates dictionary of experiment information
dlg = gui.DlgFromDict(dictionary = expInfo, title = expName)    #creates popup at beginning of experiment that asks for participant number
if dlg.OK == False: #says, if you hit escape/click cancel when that popup appears, then don't run the experiment; if this if statement didn't exist, experiment would run regardly of whether you hit escape/click cancel
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s/%s' % (expInfo['participant'], expInfo['participant'])    #creates data file name
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

flash_size = 1.15
cross_size = .07
square_mult = 528

sides_x = int(square_mult / 2)
sides_y = int(square_mult * math.sqrt(3) / 6)
bottom_y = int(-square_mult * math.sqrt(3) / 4)


##-----------------------------SHAPE COLORS-----------------------------------##

shapecolor = [.05, .05, .05]
crosscolor = [1, 1, 1]
lilcolor = [-.35, -.35, -.35]

opacity = .14


##---------------------------SOUND NOISE & LOUDNESS---------------------------##

sound_clip = sound.Sound('A')
soundfiles = [os.path.join('stimuli','ding.wav'),
              os.path.join('stimuli','chord.wav')]


##-----------------------TARGET INTERVALS & DURATION SIZES--------------------##

framelength = win.monitorFramePeriod

stagger = int(round(.0167/framelength)) #frame rate is the units
blockdelay = int(round(1.5/framelength)) #creates a slight delay after the instruction screen and before the start of each block so the onset of the block's first trial isn't too sudden
lilduration = int(round(.0333/framelength))
flash_duration = int(round(.0333/framelength))
square_start_min = int(round(1/framelength))
square_start_max = int(round(1.2/framelength))
flash_start_min = int(round(.4/framelength))
flash_start_max = int(round(.8/framelength))
lilafterflash_constant = int(round(.3/framelength))
squareafterlil = int(round(1/framelength))


##--------------------------------CREATE TIMERS-------------------------------##

globalClock = core.MonotonicClock()  # to track the time since experiment started
trialClock = core.Clock() #unlike globalclock, gets reset each trial


##---------------------------TRIAL MATRIX & RANDOMIZATION---------------------##

def round_to_multiple(f, g): #so that the num of bins can be split among even-length blocks
    return math.ceil(f * 1.0/g) * g
blocks = 1  #this line is a placeholder; change blocksreal to change # of blocks that appear
blocksreal = 8
intervals = 48
intervals = round_to_multiple(intervals, blocksreal)
reps = 12
liltrials = intervals * reps

percent_catch = .1      #if value gets too high script will quit out, since the step in the range function will equal 0 which isn't allowed

def round_to_three(x):
    return int(round(x / 3.) * 3)

catchtrials = round_to_three((percent_catch * liltrials) / (1 - percent_catch))

trials = liltrials + catchtrials
trialsperblock = trials/blocksreal

validity = .7

ptrials = 10
qtrials = 30
qcutoff = 35
startThresh = .85
acc_aim = .65


##-------------------------------INSTRUCTION SCREEN---------------------------##

inst1a = visual.TextStim(
    win = win, text = "This study will begin with " + str(ptrials + qtrials) + " practice trials. On each practice trial, like each trial in the main experiment, you will see three big squares appear on the screen.\n\nPress space to continue.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst1b = visual.TextStim(
    win = win, text = "One of these squares will then bulge shortly thereafter.\n\nPress space to continue or \"B\" to go back.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst1c = visual.TextStim(
    win = win, text = "Then on some trials you will see a little square inside one of the three squares, and on other trials you will not see any little square square. " + str(int(validity * 100)) + "% of little squares will occur at the same square as the bulge, and you are encouraged to use this information to aid your performance.\n\nPress space to continue or \"B\" to go back.", units='deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst2 = visual.TextStim(
    win = win, text = "When the little square appears, you will have one second to indicate whether it was on the left, bottom, or right of the screen by pressing, respectively, the \"A\", \"B\", or \"L\" keys. If you do not see a little square, indicate this by not pressing any button. Please gaze at the cross in the center of the screen the whole time and detect the little square with your peripheral vision.\n\nPress space to continue or \"B\" to go back.",
    units = 'deg', pos = (-8, 0), height = 1, wrapWidth = 18)

inst3 = visual.TextStim(
    win = win, text = "Again, these are practice trials, and if you do well enough on the practice trials, you will move on directly to the main experiment; otherwise you'll get a second chance to improve on the practice. As a reminder, when the little square appears, you will have one second to indicate whether it was on the left, bottom, or right of the screen by pressing, respectively, the \"A\", \"B\", or \"L\" keys.\n\nPress space to begin or \"B\" to go back.",
    units = 'deg', height = 1, wrapWidth = 20)

inst45 = visual.TextStim(
    win = win, text = "Please see the experimenter.",
    units = 'deg', height = 1, wrapWidth = 20)

inst5 = visual.TextStim(
    win = win, text = "Welcome to the beginning of the main experiment. This experiment will last about 35 minutes. It will feature trials split among " + str(blocksreal - 1) + " breaks (the breaks will be self-timed, so you can take as long as you'd like during them before proceeding to the subsequent trials).\n\nPress space to continue.",
    units = 'deg', height = 1, wrapWidth = 20)

inst6 = visual.TextStim(
    win = win, text = "The forthcoming trials will work just like the practice trials: three big squares, a bulge, and then one second to push a button to answer where the little square appeared, or not to push a button in case you did not see a little square. Again, please gaze at the cross in the center of the screen the whole time and only use peripheral vision to detect the little square.\n\nPress space to begin or \"B\" to go back.",
    units = 'deg', height = 1, wrapWidth = 20)

inst8 = visual.TextStim(
    win = win, text = "Thank you so much for your participation! Let the experimenter know that you're finished, and he'll set up the 1-minute, post-study demographic survey.",
    units = 'deg', height = 1, wrapWidth = 20)


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

outer_frame = visual.Rect(
    win = win, units = 'pix', height = square_size * 1.1, width = square_size * 1.1,
    pos = (0, 0), lineColor = shapecolor,
    fillColor = shapecolor, depth = -1.0)

lilsquare = visual.Rect(
    win = win, units = 'pix', height = lilsize, width = lilsize,
    lineColor = shapecolor, fillColor = lilcolor, depth = -1.0)

task_diagram_big_squares = visual.ImageStim(
    win = win, image = os.path.join('stimuli', 'stim_presentation_big_squares_Exp_2.png'),
    pos = (.54, 0), size = (.8, 1), texRes = 256)

task_diagram_flash = visual.ImageStim(
    win = win, image = os.path.join('stimuli', 'stim_presentation_bulge_Exp_2.png'),
    pos = (.54, 0), size = (.8, 1), texRes = 256)

task_diagram_lilsquare = visual.ImageStim(
    win = win, image = os.path.join('stimuli', 'stim_presentation_lilsquare_Exp_2.png'),
    pos = (.54, 0), size = (.8, 1), texRes = 256)

task_diagram_response = visual.ImageStim(
    win = win, image = os.path.join('stimuli', 'stim_presentation_response_Exp_2.png'),
    pos = (.54, 0), size = (.8, 1), texRes = 256)




##-----------------------CREATE EXPERIMENT MATRIX-----------------------------##

extra_valids = int(round_to_three(validity * liltrials - (6 * intervals)) / 3) # for our desired percent of valid trials, some `reps` had to include both valid and invalid trials-and then we randomized the 48 intervals assigned to these two 'mixed' reps
extra_invalids = (intervals - extra_valids) / 2

def lastreps(x, y, z):
    return np.asarray([x] * extra_valids + [y, z] * extra_invalids)
def invalids(x, y):
    return np.asarray([x, y] * (intervals / 2))

x_list = list(range(-1, 2))
y_list = [sides_y, bottom_y, sides_y]

target_x = np.concatenate([np.repeat(x_list, intervals * 2), np.repeat(x_list, intervals), np.repeat(x_list, intervals), np.repeat(x_list, int(catchtrials / 3))]) #side of screen of lil

target_y = np.concatenate([np.repeat(y_list, intervals * 2), np.repeat(y_list, intervals), np.repeat(x_list, intervals), np.repeat(y_list, int(catchtrials / 3))]) #side of screen of lil

flash_x = np.concatenate([np.repeat(x_list, intervals * 2), lastreps(-1, 0, 1), lastreps(0, -1, 1), lastreps(1, 0, -1), invalids(0, 1), invalids(-1, 1), invalids(-1, 0), np.repeat(x_list, int(catchtrials / 3))])

flash_y = np.concatenate([np.repeat(y_list, intervals * 2), lastreps(sides_y, bottom_y, sides_y), lastreps(bottom_y, sides_y, sides_y), lastreps(sides_y, bottom_y, sides_y), invalids(bottom_y, sides_y), invalids(sides_y, sides_y), invalids(sides_y, bottom_y), np.repeat(y_list, int(catchtrials / 3))])

intervals_range = list(range(0, intervals * stagger, stagger))
lil_timing = list(chain.from_iterable([random.sample(intervals_range, len(intervals_range)) for x in list(range(reps))])) # randomizes order of CTI's for each block of 48; so all 48 appear before the next rep, but the order for each 48 is random from rep to rep

if catchtrials > intervals:  # spaces out when the lilsquare comes on after the flash for catch trials
    catch_timing = [round(x * intervals * 1.0 / (catchtrials - intervals)) for x in list(range(0, (catchtrials - intervals)))]
    catch_timing += list(range(0, intervals * stagger, stagger))
else:
    catch_timing = [round(x * intervals * 1.0 / (catchtrials)) for x in list(range(0, (catchtrials)))]
np.random.shuffle(catch_timing)

expmatrix = [target_x,
            target_y,
            lil_timing + catch_timing,
            [opacity] * liltrials + [0] * catchtrials, #opacity of lil
            flash_x, # side of screen of flash
            flash_y, # y-value of flash
            np.asarray([random.randrange(sides_x * target_x[i] - square_size / 2 + lilsize / 2, sides_x * target_x[i] + square_size / 2 - lilsize / 2) for i in list(range(trials))]), # randomize target's x value
            np.asarray([random.randrange(target_y[i] - square_size / 2 + lilsize / 2, target_y[i] + square_size / 2 - lilsize / 2) for i in list(range(trials))])] # randomize target's y value
print expmatrix
print str(len(expmatrix[0])) + str(len(expmatrix[1])) + str(len(expmatrix[2])) + str(len(expmatrix[3])) + str(len(expmatrix[4])) + str(len(expmatrix[5])) + str(len(expmatrix[6])) + str(len(expmatrix[7]))

#randomization sequence
randomseq = list(range(int(trials)))
np.random.shuffle(randomseq)







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
        inst1b.setAutoDraw(False)
        inst1c.setAutoDraw(False)
        inst2.setAutoDraw(False)
        inst3.setAutoDraw(False)
        task_diagram_flash.setAutoDraw(False)
        task_diagram_lilsquare.setAutoDraw(False)
        task_diagram_response.setAutoDraw(False)
        inst1a.setAutoDraw(True)
        task_diagram_big_squares.setAutoDraw(True)
    elif advance == 1:
        inst1a.setAutoDraw(False)
        inst1c.setAutoDraw(False)
        inst2.setAutoDraw(False)
        inst3.setAutoDraw(False)
        task_diagram_big_squares.setAutoDraw(False)
        task_diagram_lilsquare.setAutoDraw(False)
        task_diagram_response.setAutoDraw(False)
        inst1b.setAutoDraw(True)
        task_diagram_flash.setAutoDraw(True)
    elif advance == 2:
        inst1b.setAutoDraw(False)
        inst2.setAutoDraw(False)
        inst3.setAutoDraw(False)
        task_diagram_flash.setAutoDraw(False)
        task_diagram_response.setAutoDraw(False)
        inst1c.setAutoDraw(True)
        task_diagram_lilsquare.setAutoDraw(True)
    elif advance == 3:
        inst1c.setAutoDraw(False)
        inst3.setAutoDraw(False)
        task_diagram_lilsquare.setAutoDraw(False)
        inst2.setAutoDraw(True)
        task_diagram_response.setAutoDraw(True)
    elif advance == 4:
        inst2.setAutoDraw(False)
        task_diagram_response.setAutoDraw(False)
        inst3.setAutoDraw(True)
    else:
        inst3.setAutoDraw(False)


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
        np.random.shuffle(randomseq)  # reshuffle the order of trials so that practice/staircase trials are not in the same order as experimental trials
    while q_acc < qcutoff:
        acclist = []
        if rep == 1:


            ##------TELL PTCPT TO SEE EXPERIMENTER BEFORE 2ND RESTART---------##

            while qloop == 2:
                inst45.setAutoDraw(True)
                if(event.getKeys(keyList = ["space"])):
                    inst45.setAutoDraw(False)
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
                            inst6.setAutoDraw(False)
                            inst5.setAutoDraw(True)
                        elif (advance == 1):
                            inst5.setAutoDraw(False)
                            inst6.setAutoDraw(True)
                        else:
                            inst6.setAutoDraw(False)
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

                cross.setAutoDraw(True) # cross begins on screen; is located outside of while loop since it is on screen the entire trial
                frameN = -1
                continueRoutine = True
                lenkeylist = 0
                overalltime = globalClock.getTime()
                routineTimer = core.Clock()


                ##----------------------SET TARGET OPACITY--------------------##

                if rep == 0:
                    trialopacity = opacity * startThresh
                elif rep == 1:
                    trialopacity = opacity * trial # trial here refers to the updated value in the staircase, not the actual trial number
                    trial = ptrials
                    ptrials += 1
                else:
                    extracted_opacity = expmatrix[3][randomseq[trial]]
                    if extracted_opacity > 0:
                        if noncatch_count > 0 and noncatch_count % 16 == 0:
                            repstair_avg = sum(repstaircase) * 1.0 / len(repstaircase)
                            q_opacity += .033 * (acc_aim - repstair_avg)
                            repstaircase = []
                        noncatch_count += 1
                    trialopacity = extracted_opacity * q_opacity #whether no opacity * threshold

                lilsquare.setOpacity(trialopacity)


                ##----------------SET TARGET & FLASH LOCATIONS----------------##

                lil_side = expmatrix[0][randomseq[trial]]
                flash_side = expmatrix[4][randomseq[trial]]
                outer_frame.setPos([flash_side * sides_x, expmatrix[5][randomseq[trial]]])
                lilsquare.setPos([[expmatrix[6][randomseq[trial]], expmatrix[7][randomseq[trial]]]])


                ##--------------SET START & DURATION OF STIMULI---------------##

                square_start = random.randrange(square_start_min, square_start_max + 1)
                flash_start = square_start + random.randrange(flash_start_min, flash_start_max + 1)
                flash_end = flash_start + flash_duration

                lilafterflash = expmatrix[2][randomseq[trial]]
                lilstart = flash_end + lilafterflash_constant + lilafterflash
                lilend = lilstart + lilduration

                squareend = lilend + squareafterlil
                square_duration = squareend - square_start
                trial_duration = square_start + square_duration


                ##-------------------RESET TRIAL CLOCK------------------------##

                trialClock.reset()  # clock



                ##------------------------------------------------------------##
                ##------------------WHILE LOOP BEGINS-------------------------##
                ##------------------------------------------------------------##


                while continueRoutine and frameN <= trial_duration:
                    if event.getKeys(keyList = ["escape"]):
                        core.quit()
                    # get current time
                    t = trialClock.getTime()
                    frameN += 1  # number of completed frames (so 0 is the first frame)

                    key = event.getKeys(keyList = ["l","L","a","A", "b", "b"])

                    ##----------------SHAPES UPDATE---------------------------##

                    if frameN == square_start:
                        square_right.tStart = t
                        square_right.frameStart = frameN
                        square_right.setAutoDraw(True)
                        square_left.setAutoDraw(True)
                        square_bottom.setAutoDraw(True)
                    elif frameN == flash_start:
                        outer_frame.tStart = t
                        outer_frame.frameStart = frameN
                        outer_frame.setAutoDraw(True)
                    elif frameN == flash_end:
                        outer_frame.setAutoDraw(False)
                        outer_frame.tEnd = t
                        outer_frame.frameEnd = frameN
                    elif frameN == lilstart:
                        lilsquare.tStart = t
                        lilsquare.frameStart = frameN
                        lilsquare.setAutoDraw(True)
                    elif frameN == lilend:
                        lilsquare.setAutoDraw(False)
                        lilsquare.tEnd = t
                        lilsquare.frameEnd = frameN
                    elif frameN == trial_duration:
                         key = ['nope']
                    if len(key) > 0:
                        if frameN < lilstart:
                            lenkeylist += 1
                        else:
                            continueRoutine = False
                            square_right.setAutoDraw(False)
                            square_left.setAutoDraw(False)
                            square_bottom.setAutoDraw(False)
                            cross.setAutoDraw(False)


                            ##-------------CHECK FOR RESPONSE-----------------##

                            if (key == ['nope'] and trialopacity == 0) or ((key == ['l'] and lil_side == 1) or (key == ['a'] and lil_side == -1) or (key == ['b'] and lil_side == 0)):
                                acc = 1
                                soundp = soundfiles[0]
                            else:
                                acc = 0
                                soundp = soundfiles[1]
                            if rep == 1:
                                trials.addResponse(acc)
                            acclist.append(acc)


                    event.clearEvents()     # Clear the previously pressed keys; too-early key presses will automatically register as incorrect


                    ##--------CHECK ALL IF COMPONENTS HAVE FINISHED-----------##

                    if not continueRoutine:  # a component has requested a forced-end of Routine
                        break
                    else:
                        win.flip()


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
                thisExp.addData('flash_circleStartTime', outer_frame.tStart)
                thisExp.addData('flash_circleStartFrame', outer_frame.frameStart)
                thisExp.addData('flash_circleEndTime', outer_frame.tEnd)
                thisExp.addData('flash_circleEndFrame', outer_frame.frameEnd)
                thisExp.addData('lilsquareStartTime', lilsquare.tStart)
                thisExp.addData('lilsquareStartFrame', lilsquare.frameStart)
                thisExp.addData('lilsquareEndTime', lilsquare.tEnd)
                thisExp.addData('lilsquareEndFrame', lilsquare.frameEnd)
                thisExp.addData('lilafterflash', lilafterflash)
                thisExp.addData('ButtonPressTime', t)
                thisExp.addData('Key', key[0])
                thisExp.addData('false_presses', lenkeylist)
                thisExp.addData('Acc', acc)
                thisExp.addData('Opacity', trialopacity)
                thisExp.addData('FlashSide', flash_side)
                thisExp.addData('CorrSide', lil_side)
                thisExp.nextEntry()
                sound_clip.setSound(soundp)
                sound_clip.play() #correct and incorrect sounds are both 250 ms

            if rep == 1:
                q_opacity = trials.mean()
                q_acc = round(np.mean(acclist) * 100)
                qloop += 1



inst8.setAutoDraw(True)
while len(event.getKeys(keyList=['space'])) == 0:
    win.flip()
inst8.setAutoDraw(False)

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename + '.csv')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
