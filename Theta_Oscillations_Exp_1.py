# PsychoPy is required for this experiment
from psychopy import locale_setup, sound, gui, visual, core, data, event, logging, clock
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy.random import random, randint, normal, shuffle
import random
import os  # handy system and path functions
import sys  # to get file system encoding
import math

# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__)).decode(sys.getfilesystemencoding())
os.chdir(_thisDir)

# Store info about the experiment session
expName = 'Theta_Oscillations_Exp_1.py'
expInfo = {'participant': ''}     #creates dictionary of experiment information
dlg = gui.DlgFromDict(dictionary=expInfo, title=expName)    #creates popup at beginning of experiment that asks for participant number
if dlg.OK == False: #says, if you hit escape/click cancel when that popup appears, then don't run the experiment; if this if statement didn't exist, experiment would run regardly of whether you hit escape/click cancel
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s' % (expInfo['participant'], expName)    #creates data file name
# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(extraInfo=expInfo, dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file; unclear what this line does




##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##
##-----------------------------START CODE-------------------------------------##
##----------------------------------------------------------------------------##
##----------------------------------------------------------------------------##



# Setup the Window
win = visual.Window(
    size=(1024, 768), color=[.15,.15,.15],fullscr=True, allowGUI=False,
    monitor='testMonitor', useFBO=True)
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

circle_radius = [12,12]
cross_size = .07
square_x = 280
square_y = 0
flash_dist = 130


##-----------------------------SHAPE COLORS-----------------------------------##

shapecolor = [.05,.05,.05]
crosscolor = [1,1,1]
flashcolor = [.52,.52,.52]
lilcolor = [-.35,-.35,-.35]

opacity = .14


##-----------------------TARGET INTERVALS & DURATION SIZES--------------------##

framelength = win.monitorFramePeriod

stagger = int(round(.0167/framelength)) #frame rate is the units
blockdelay = int(round(1.5/framelength)) #creates a slight delay after the instruction screen and before the start of each block so the onset of the block's first trial isn't too sudden
lilduration = int(round(.0333/framelength))
flash_duration = int(round(.0333/framelength))
square_start_min = int(round(1/framelength))
square_start_max = int(round(1.2/framelength))
flash_start_min = int(round(.4167/framelength))
flash_start_max = int(round(.8333/framelength))
lilafterflash_constant = int(round(.5/framelength))
squareafterlil = int(round(1/framelength))


##--------------------------------CREATE TIMERS-------------------------------##

globalClock = core.MonotonicClock()  # to track the time since experiment started
trialClock = core.Clock() #unlike globalclock, gets reset each trial


##---------------------------TRIAL MATRIX & RANDOMIZATION---------------------##

def round_to_multiple(f,g): #so that the num of bins can be split among even-length blocks
    return math.ceil(f * 1.0 /g) * g
blocks = 1  #this line is a placeholder; change blocksreal to change # of blocks that appear
blocksreal = 8
intervals = 80
intervals = round_to_multiple(intervals, blocksreal)
reps = 4
liltrials = intervals * reps

catch_to_noncatch = 10      #if value gets too high script will quit out, since the step in the range function will equal 0 which isn't allowed
catchtrials = int(liltrials*catch_to_noncatch/100.0)

trials = liltrials + catchtrials
trialsperblock = trials/blocksreal

ptrials = 10
qtrials = 30
qcutoff = 35
startThresh = .85


##-------------------------------INSTRUCTION SCREEN---------------------------##

inst1a = visual.TextStim(
    win=win, text = "This study will begin with " + str(ptrials + qtrials) + " practice trials. On each practice trial, like each trial in the main experiment, you will see two big squares appear on the screen.\n\nPress space to continue.",
    units='deg', pos=(-8, 0), height = 1, wrapWidth = 18)

inst1b = visual.TextStim(
    win=win, text = "You will then see four white circles flash shortly thereafter, and you can ignore this flash as it has no behavioral relevance to the task.\n\nPress space to continue or \"B\" to go back.",
    units='deg', pos=(-8, 0), height = 1, wrapWidth = 18)

inst1c = visual.TextStim(
    win=win, text = "Then on some trials you will see a little square inside one of the two squares, and on other trials you will not see any little square.\n\nPress space to continue or \"B\" to go back.", units='deg', pos=(-8, 0), height = 1, wrapWidth = 18)

inst2 = visual.TextStim(
    win=win, text = "When the little square appears, you will have one second to indicate whether it was on the left or right of the screen by pressing, respectively, the \"A\" or \"L\" keys. If there is no little square on the screen, you do not have to press any key. Please gaze at the cross in the center of the screen the whole time and detect the little square with your peripheral vision. Finally, please stay in the head mount during throughout the experiment.\n\nPress space to continue or \"B\" to go back.",
    units='deg',  pos=(-8, 0), height = 1,wrapWidth = 18)

inst3 = visual.TextStim(
    win=win, text="Again, these are practice trials, and if you get more than " + str(qcutoff) + "% of the trials correct you will move on directly to the main experiment; otherwise you'll get a second chance to improve on the practice.\n\nPress space to begin or \"B\" to go back.",
    units='deg', height = 1, wrapWidth = 20)

inst45 = visual.TextStim(
    win=win, text="Please see the experimenter.",
    units='deg', height = 1, wrapWidth = 20)

inst5 = visual.TextStim(
    win=win, text="Welcome to the beginning of the main experiment. This experiment will last about 20 minutes. It will feature " + str(trials) + " trials split among " + str(blocksreal - 1) + " breaks (the breaks will be self-timed, so you can take as long as you'd like during them before proceeding to the subsequent trials).\n\nPress space to continue.",
    units='deg', height = 1, wrapWidth = 20)

inst6 = visual.TextStim(
    win=win, text="The forthcoming trials will work just like the practice trials: two big squares, a flash, and then one second to push a button to answer on which side the little square appeared, or not to push a button in case you did not see a little square. Again, please gaze at the cross in the center of the screen the whole time and only use peripheral vision to detect the little square.\n\nPress space to begin or \"B\" to go back.",
    units='deg', height = 1, wrapWidth = 20)

inst8 = visual.TextStim(
    win=win, text="Thank you so much for your participation! Let the experimenter know that you're finished, and he'll set up the 1-minute, post-study demographic survey.",
    units='deg', height = 1, wrapWidth = 20)


##----------------------MAKE THE SEARCH ARRAY ITEMS---------------------------##

cross = visual.ShapeStim(
    win=win, vertices=[[0,0],[0, .2306 * cross_size], [0,0],
    [0, -.2306 * cross_size],[0,0],[.1501 * cross_size, 0],
    [0,0],[-.1501 * cross_size,0],[0,0]], lineWidth=2,
    lineColor=crosscolor, depth=-1.0)

square_right = visual.Rect(
    win=win, units = 'pix', height = square_size,
    width = square_size, pos=(square_x,square_y),
    lineColor=shapecolor, fillColor=shapecolor,
    depth=-1.0)

square_left = visual.Rect(
    win=win, units = 'pix', height = square_size,
    width = square_size, pos=(-square_x,square_y),
    lineColor=shapecolor, fillColor=shapecolor,
    depth=-1.0)

flash_circle_outer = visual.Circle(
    win=win, radius=circle_radius, units="pix",
    edges=2048, fillColor=flashcolor,
    lineColor=flashcolor)

flash_circle_inner = visual.Circle(
    win=win, radius=circle_radius, units="pix",
    edges=2048, fillColor=flashcolor,
    lineColor=flashcolor)

flash_circle_top = visual.Circle(
    win=win, radius=circle_radius, units="pix",
    edges=2048, fillColor=flashcolor,
    lineColor=flashcolor)

flash_circle_bottom = visual.Circle(
    win=win, radius=circle_radius, units="pix",
    edges=2048, fillColor=flashcolor,
    lineColor=flashcolor)

lilsquare = visual.Rect(
    win=win, units = 'pix', height = lilsize,
    width = lilsize, lineColor=shapecolor,
    fillColor=lilcolor, depth=-1.0)

task_diagram_big_squares = visual.ImageStim(
    win=win, image= os.path.join('stimuli','stim_presentation_big_squares.png'),
    pos=(.54, 0), size=(.8, 1),
    texRes=256)

task_diagram_flash = visual.ImageStim(
    win=win, image= os.path.join('stimuli', 'stim_presentation_flash.png'),
    pos=(.54, 0), size=(.8, 1),
    texRes=256)

task_diagram_lilsquare = visual.ImageStim(
    win=win, image= os.path.join('stimuli', 'stim_presentation_lilsquare.png'),
    pos=(.54, 0), size=(.8, 1),
    texRes=256)

task_diagram_response = visual.ImageStim(
    win=win, image= os.path.join('stimuli', 'stim_presentation_response.png'),
    pos=(.54, 0), size=(.8, 1),
    texRes=256)


##-------------------------JITTER THE STAIRCASE-------------------------------##

jitterlevel1 = .05
jitterlevel2 = .1
jitterlevel1_freq = .2 #meaning .2 freq for level1 above threshold, and another .2 freq for level1 below threshold
jitterlevel2_freq = .05 #ditto as comment for the line above

def jitter():
    return (np.repeat(np.arange(1 - jitterlevel1, 1 + jitterlevel1*2, jitterlevel1*2),
    liltrials*jitterlevel1_freq), np.repeat(np.arange(1 - jitterlevel2, 1 + jitterlevel2*2, jitterlevel2*2), liltrials*jitterlevel2_freq))

jittered_array = np.concatenate([jitter()[0], jitter()[1], np.repeat(1, liltrials - (len(jitter()[0]) + len(jitter()[1])))])  #trial without jitter is just the number of trials - the number of trials that will be jittered
np.random.shuffle(jittered_array) #expmatrix never gets shuffled, just randomseq; this line is needed so random relation between jitter and expmatrix


##-----------------------CREATE EXPERIMENT MATRIX-----------------------------##

expmatrix = [np.concatenate([np.repeat(range(-1,2,2),int(liltrials/2)), np.repeat(range(-1,2,2),int(catchtrials/2))]), #side of screen of lil
                reps * range(0,intervals * stagger,stagger) + [round(x * intervals*1.0/catchtrials) + 1 for x in range(0, catchtrials)], #part after the + spaces out when the lilsquare comes on after the flash for catch trials
                [opacity]*liltrials + [0]*catchtrials, #opacity of lil
                np.concatenate([np.repeat([-1,1,-1,1],int(liltrials/4)), np.repeat([-1,1],int(catchtrials/2))]), # side of screen of flash
                np.concatenate((jittered_array, np.repeat(0, catchtrials)))]

#randomization sequence
randomseq = range(int(trials))
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
    if event.getKeys(keyList=["space"]):
        advance += 1
    elif event.getKeys(keyList=["b"]):
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



trials = range(ptrials)
q_opacity = 0
qloop = 0

for rep in range(3):
    q_acc = 0
    if rep == 2:
        blocks = blocksreal
        np.random.shuffle(randomseq)
    while q_acc < qcutoff:
        acclist = []
        if rep == 1:


            ##------TELL PTCPT TO SEE EXPERIMENTER BEFORE 2ND RESTART---------##

            while qloop == 2:
                inst45.setAutoDraw(True)
                if(event.getKeys(keyList=["space"])):
                    inst45.setAutoDraw(False)
                    qloop = 1
                win.flip()


            ##---SHOW INSTRUCTIONS AGAIN IF PTCPT HASN'T REACHED THRESHOLD----##

            if qloop == 1:
                inst4 = visual.TextStim(
                    win=win, text=str(q_acc) + "% of your responses were correct, which is less than the " + str(qcutoff) + "% threshold. As a reminder, when the little square appears, you will have one second to indicate whether it was on the left or right of the screen by pressing, respectively, the \"A\" or \"L\" keys. If there is no little square on the screen, you do not have to press any key. Finally, please detect the little square with your peripheral vision and gaze at the cross in the center of the screen the whole time.\n\nPlease press the space bar to try again.",
                    units='deg', height = 1, wrapWidth = 20)
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
                pThreshold=.62, gamma=0.05, #.82, says http://www.psychopy.org/api/data.html#psychopy.data.QuestHandler is "equivalent to a 3 up 1 down standard staircase"
                nTrials=qtrials, minVal= .01, maxVal= 4)
        else:
            q_acc = qcutoff + 1
            frameN = -1 # number of completed frames (so 0 is the first frame)
            while frameN < blockdelay: # creates a 1.5 second blank screen that appears after the instructions but before the first practice trial begins
                frameN += 1
                win.flip()

        for block in range(blocks):
            if rep == 2:
                continueRoutineInst = True

                inst7 = visual.TextStim(
                    win=win, text= "You've reached break " + str(block) + " of " + str(blocks-1) + ". This break is self-timed, so whenever you're ready press spacebar to continue the study. Feel free to leave the head mount during the break, but when the break ends please return to your original position on the head mount.",
                    units='deg', height = 1, wrapWidth = 20)

                advance = 0 # a variable that advances the instruction screen, as well as lets them go back to see a previous instruction screen

                while continueRoutineInst:
                    if block == 0:
                        if(event.getKeys(keyList=["space"])):
                            advance += 1
                        elif(event.getKeys(keyList=["b"])):
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
                        if (event.getKeys(keyList=["space"])):
                            continueRoutineInst = False
                            inst7.setAutoDraw(False)

                    if continueRoutineInst:
                        win.flip()

                frameN = -1
                while frameN < blockdelay:
                    frameN += 1
                    win.flip()

                startingtrial = block*trialsperblock
                trials = range(startingtrial,startingtrial + trialsperblock)



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
                    trialopacity = opacity*startThresh
                elif rep == 1:
                    trialopacity = opacity*trial
                    trial = ptrials
                    ptrials += 1
                else:
                    trialopacity = expmatrix[2][randomseq[trial]]*q_opacity*expmatrix[4][randomseq[trial]] #whether no opacity * threshold * potentially jittering

                lilsquare.setOpacity(trialopacity)


                ##----------------SET TARGET & FLASH LOCATIONS----------------##

                flash_side = expmatrix[3][randomseq[trial]]
                flash_circle_outer.setPos([(square_x + flash_dist) * flash_side, square_y])
                flash_circle_inner.setPos([(square_x - flash_dist) * flash_side, square_y])
                flash_circle_top.setPos([square_x * flash_side, square_y + flash_dist])
                flash_circle_bottom.setPos([square_x * flash_side,square_y - flash_dist])
                lil_x = random.randrange(square_x - square_size/2 + lilsize/2, square_x + square_size/2 - lilsize/2)
                lil_side = expmatrix[0][randomseq[trial]]
                lil_x = lil_x * lil_side
                lil_y = random.randrange(square_y - square_size/2 + lilsize/2, square_y + square_size/2 - lilsize/2)
                lilsquare.setPos([lil_x, lil_y])


                ##--------------SET START & DURATION OF STIMULI---------------##

                square_start = random.randrange(square_start_min, square_start_max + 1)
                flash_start = square_start + random.randrange(flash_start_min, flash_start_max + 1)
                flash_end = flash_start + flash_duration

                lilafterflash = expmatrix[1][randomseq[trial]]
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
                    if event.getKeys(keyList=["escape"]):
                        core.quit()
                    # get current time
                    t = trialClock.getTime()
                    frameN += 1  # number of completed frames (so 0 is the first frame)

                    key = event.getKeys(keyList=["l","L","a","A"])

                    ##----------------SHAPES UPDATE---------------------------##

                    if frameN == square_start:
                        square_right.tStart = t
                        square_right.frameStart = frameN
                        square_right.setAutoDraw(True)
                        square_left.setAutoDraw(True)
                    elif frameN == flash_start:
                        flash_circle_outer.tStart = t
                        flash_circle_outer.frameStart = frameN
                        flash_circle_outer.setAutoDraw(True)
                        flash_circle_inner.setAutoDraw(True)
                        flash_circle_top.setAutoDraw(True)
                        flash_circle_bottom.setAutoDraw(True)
                    elif frameN == flash_end:
                        flash_circle_outer.setAutoDraw(False)
                        flash_circle_inner.setAutoDraw(False)
                        flash_circle_top.setAutoDraw(False)
                        flash_circle_bottom.setAutoDraw(False)
                        flash_circle_outer.tEnd = t
                        flash_circle_outer.frameEnd = frameN
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
                            cross.setAutoDraw(False)


                            ##-------------CHECK FOR RESPONSE-----------------##

                            if (key == ['nope'] and trialopacity == 0) or ((key == ['l'] and lil_side == 1) or (key == ['a'] and lil_side == -1)):
                                acc = 1
                            else:
                                acc = 0
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
                    thisExp.addData('Trial', trial)
                thisExp.addData('ButtonPressTimeinOverallExp', overalltime)
                thisExp.addData('squareStartTime', square_right.tStart)
                thisExp.addData('squareStartFrame', square_right.frameStart)
                thisExp.addData('flash_circleStartTime', flash_circle_outer.tStart)
                thisExp.addData('flash_circleStartFrame', flash_circle_outer.frameStart)
                thisExp.addData('flash_circleEndTime', flash_circle_outer.tEnd)
                thisExp.addData('flash_circleEndFrame', flash_circle_outer.frameEnd)
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

            if rep == 1:
                q_opacity = trials.mean()
                q_acc = round(np.mean(acclist)*100)
                qloop += 1



inst8.setAutoDraw(True)
while len(event.getKeys(keyList=['space'])) == 0:
    win.flip()
inst8.setAutoDraw(False)

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
