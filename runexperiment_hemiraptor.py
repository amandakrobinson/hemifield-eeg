#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 14:10:15 2019

@author: tgro5258
"""

from psychopy import core, event, visual, parallel, gui
import random,sys,json,requests,os,itertools
from glob import glob
import pandas as pd
import numpy as np

# debug things
debug_testsubject = 1
debug_usedummytriggers = 1
debug_windowedmode = 1
debug_onesequence = 1


stimsize = 128 #size of stim
wordsize = 40 #height of word
stimshift = 128 #shift away from fixation

objectstimuli = sorted(glob('Stimuli/*.png'))
wordstimuli = ['fish','bird','face','boat','tree','tool']

if debug_testsubject:
    subjectnr = 0
else:
    # Get subject info
    subject_info = {'Subject number (update participants.tsv)':''}
    if not gui.DlgFromDict(subject_info,title='Enter subject info (update participants.tsv):').OK:
        print('User hit cancel at subject information')
        exit()
    try:
        subjectnr = int(subject_info['Subject number (update participants.tsv)'])
    except:
        raise

outfn = 'sub-%02i_task-rsvp_events.csv'%subjectnr
if not debug_testsubject and os.path.exists(outfn):
    raise Exception('%s exists'%outfn)

nstimuli = 48

nsequence = 1 if debug_onesequence else 12

random.seed(subjectnr)

refreshrate = 60
feedbackduration = .5 - .5/refreshrate
fixationduration = 1 - .5/refreshrate
stimduration = .1 - .5/refreshrate
isiduration = .2 - .5/refreshrate

trigger_stimon = 1
trigger_stimoff = 2
trigger_sequencestart = 3
trigger_duration = 0.005
trigger_port = 0xcff8 #badham
# trigger_port = 0xd050 #GT322

webhook_url='https://hooks.slack.com/services/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'

stimnum = list(range(nstimuli))

# conditions: central, (L, R), LR

# Center condition: 4 repeats per stream
eventlist1 = pd.DataFrame()
for i in range(nsequence):
    stream = []
    for j in range(4):
        x=[]
        while not x or stream and len(set(stream[-4:]).intersection(set(x[:4]))):
            x = random.sample(range(nstimuli),nstimuli)
        stream+=x
    for (x1,x2) in enumerate(stream):
        eventlist1=eventlist1.append(pd.DataFrame([[i,x1,x2]], columns=['sequencenumber','presentationnumber','stimC']), ignore_index=1)
    print('Condition 1: stream %i/%i'%(i+1,nsequence))
eventlist1['stimL']=-1
eventlist1['stimR']=-1
eventlist1['condition']=1
eventlist1['conditionlabel']='center'

# LR Simultaneous condition = all combinations of stims
eventlist2 = pd.DataFrame()
condition2list = list(itertools.product(stimnum,stimnum))
#shuffle list while we have repeats left or right
random.shuffle(condition2list)
#cut up into sequences:
for i in range(nsequence):
    stream = [condition2list.pop() for x in range(192)]
    while not np.all(np.all(np.diff(stream,axis=0),axis=0)):
        random.shuffle(stream)
    for (x1,(x2,x3)) in enumerate(stream):
        eventlist2=eventlist2.append(pd.DataFrame([[i,x1,x2,x3]], columns=['sequencenumber','presentationnumber','stimL','stimR']), ignore_index=1)
    print('Condition 2: stream %i/%i'%(i+1,nsequence))
eventlist2['stimC']=-1
eventlist2['condition']=3
eventlist2['conditionlabel']='left+right'

# left-right Condition
eventlist3 = pd.DataFrame()
eventlist4 = pd.DataFrame()
for i in range(nsequence):
    streamtotal=[]
    loctotal=[]
    for k in range(2):
        #pseudo-randomise left/right
        loc=[]
        for x in range(48):
            loc+=random.sample([0, 1],2)
        loctotal+=loc
        stream=[0,0]
        while not all(np.diff(stream)):
            streamL = []
            for j in range(2):
                x=[]
                while not x or streamL and len(set(streamL[-4:]).intersection(set(x[:4]))):
                    x = random.sample(range(nstimuli),nstimuli)
                streamL+=x
            streamR = streamL[::-1]
            stream=[]
            for (j,e) in enumerate(loc):
                if e:
                    stream.append(streamL.pop())
                else:
                    stream.append(streamR.pop())
        streamtotal+=stream
        loctotal+=loc
    for (x1,(x2,x3)) in enumerate(zip(loctotal,streamtotal)):
        eventlist3=eventlist3.append(pd.DataFrame([[i,x1,x3 if x2 else -1,x3 if not x2 else -1]], columns=['sequencenumber','presentationnumber','stimL','stimR']), ignore_index=1)
        eventlist4=eventlist4.append(pd.DataFrame([[i,x1,x3 if not x2 else -1,x3 if x2 else -1]], columns=['sequencenumber','presentationnumber','stimL','stimR']), ignore_index=1)
    print('Condition 3: stream %i/%i'%(i+1,nsequence))
eventlist3['stimC']=-1
eventlist3['condition']=2
eventlist3['conditionlabel']='left-right'
eventlist4['stimC']=-1
eventlist4['condition']=2
eventlist4['conditionlabel']='left-right'

#stack
eventlist=pd.DataFrame()
conditiontriples = []
seqorders = [random.sample(range(nsequence),nsequence) for i in range(4)]
C = [eventlist1, eventlist2, eventlist3, eventlist4]
n = -1
for x in range(nsequence):
    conditiontriple=random.sample(range(4),4)
    for c in conditiontriple:
        seq = seqorders[c][x]
        idx = C[c].sequencenumber==seq
        e = C[c][idx]
        n+=1
        e['streamnumber'] = n
        #make targets
        ntargets = random.randint(3,6)
        targetpos=[1, 1]
        t = [0 for x in range(len(e))]
        while len(targetpos)>1 and any(np.diff(targetpos)<10):
            targetpos = sorted(random.sample(range(10,len(e)-10),ntargets))
        for p in targetpos:
            t[p]=random.sample('LRC',1)[0]
        e['istarget'] = t
        
        eventlist=eventlist.append(e, ignore_index=1)

    
def writeout(eventlist):
    with open(outfn,'w') as out:
        eventlist.to_csv(out,index_label='eventnumber')

writeout(eventlist)

# =============================================================================
# %% START
# =============================================================================
try:
    if debug_windowedmode:
        win=visual.Window([700,700],units='pix',color=(-.1, -.1, -.1))
    else:
        win=visual.Window(units='pix',fullscr=True)
    mouse = event.Mouse(visible=False)

    fixation = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=15,
        name='fixation', autoLog=False, color='black')
    fixationL = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=fixation.size,
        name='fixation', autoLog=False, color='black', pos=(-stimshift,0))
    fixationR = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=fixation.size,
        name='fixation', autoLog=False, color='black', pos=(stimshift,0))
    feedback = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=fixation.size,
        name='feedback', autoLog=False)
    progresstext = visual.TextStim(win,text='',pos=(0,100),name='progresstext')
    sequencestarttext = visual.TextStim(win,text='Press button 1 or 4 to start the sequence\nPress button 2 or 3 when you see a red dot',pos=(0,50),name='sequencestarttext')

    filesep='/'
    if sys.platform == 'win32':
        filesep='\\'
        
    stimtex=[]
    for (i,y) in enumerate(objectstimuli):
        stimtex.append(visual.ImageStim(win,y,size=stimsize,name=y.split(filesep)[1]))
    for (i,y) in enumerate(wordstimuli):
        stimtex.append(visual.TextStim(win,text=y.lower(),height=wordsize,italic=0,name='%s_%i_lower'%(y,i+1)))
        stimtex.append(visual.TextStim(win,text=y.lower(),height=wordsize,italic=1,name='%s_%i_lower_it'%(y,i+2)))
        stimtex.append(visual.TextStim(win,text=y.upper(),height=wordsize,italic=0,name='%s_%i_upper'%(y,i+3)))
        stimtex.append(visual.TextStim(win,text=y.upper(),height=wordsize,italic=1,name='%s_%i_upper_it'%(y,i+4)))
        
    targetstimtex = visual.GratingStim(win, tex=None, mask='gauss', sf=0, size=fixation.size,
        name='target', autoLog=False, color='red')
    
    def send_dummy_trigger(trigger_value):
        core.wait(trigger_duration)
            
    def send_real_trigger(trigger_value):
        trigger_port.setData(trigger_value)
        core.wait(trigger_duration)
        trigger_port.setData(0)
    
    if debug_usedummytriggers:
        sendtrigger = send_dummy_trigger
    else:
        trigger_port = parallel.ParallelPort(address=trigger_port)
        trigger_port.setData(0)
        sendtrigger = send_real_trigger

    nevents = len(eventlist)
    nsequences = eventlist['streamnumber'].iloc[-1]+1
    sequencenumber = -1
    for eventnr in range(len(eventlist)):
        first = eventlist['streamnumber'].iloc[eventnr]>sequencenumber
        if first:
            
            writeout(eventlist)
            sequencenumber = eventlist['streamnumber'].iloc[eventnr]
            condition = eventlist['condition'].iloc[eventnr]
            last_target = -99
            correct=0
            
            if not debug_testsubject:
                try:
                    slack_data={'text':'sub-%02i seq %i/%i (hemiraptor) <@tijlgrootswagers> <@amanda> <@U9C24ECQ7>'%(subjectnr,sequencenumber+1,nsequences),'channel':'#eeglab','username':'python'}
                    response = requests.post(webhook_url, data=json.dumps(slack_data),headers={'Content-Type': 'application/json'})
                except:
                    pass

            progresstext.text = '%i / %i'%(1+sequencenumber,nsequences)
            progresstext.draw()
            sequencestarttext.draw()
            fixation.draw()
            fixationL.draw()
            fixationR.draw()
            win.flip()
            k=event.waitKeys(keyList='afq', modifiers=False, timeStamped=True)
            if k[0][0]=='q':
                raise Exception('User pressed q')
            fixation.draw()
            fixationL.draw()
            fixationR.draw()
            time_fixon = win.flip()
            sendtrigger(trigger_sequencestart)
            while core.getTime() < time_fixon + fixationduration:pass
        
        response=0
        rt=0
        fixation.draw()
        fixationL.draw()
        fixationR.draw()
        stimnameL,stimnameR,stimnameC='','',''
        if eventlist['stimL'].iloc[eventnr]>-1:
            stim = stimtex[eventlist['stimL'].iloc[eventnr]]
            stim.pos = (-stimshift,0)
            stim.draw()
            stimnameL = stim.name
        if eventlist['stimR'].iloc[eventnr]>-1:
            stim = stimtex[eventlist['stimR'].iloc[eventnr]]
            stim.pos = (stimshift,0)
            stim.draw()
            stimnameR = stim.name
        if eventlist['stimC'].iloc[eventnr]>-1:
            stim = stimtex[eventlist['stimC'].iloc[eventnr]]
            stim.pos = (0,0)
            stim.draw()
            stimnameC = stim.name
            
        time_stimon=win.flip()
        sendtrigger(trigger_stimon)
        while core.getTime() < time_stimon + stimduration:pass
        fixation.draw()
        fixationL.draw()
        fixationR.draw()
        if eventlist['istarget'].iloc[eventnr]:
            if eventlist['istarget'].iloc[eventnr]=='L':
                targetstimtex.pos=(-stimshift,0)
            if eventlist['istarget'].iloc[eventnr]=='R':
                targetstimtex.pos=(stimshift,0)
            if eventlist['istarget'].iloc[eventnr]=='C':
                targetstimtex.pos=(0,0)
            targetstimtex.draw()
        time_stimoff=win.flip()
        sendtrigger(trigger_stimoff)
        
        if eventlist['istarget'].iloc[eventnr]:
            last_target=time_stimoff
        
        correct=0
        k=event.getKeys(keyList='sdq', modifiers=False, timeStamped=True)
        if k:
            response=k[0][0]
            rt=k[0][1]
            if response=='q':
                raise Exception('User pressed q')
            else:
                response=1
            correct = rt-last_target < 1
        
        eventlist.at[eventnr, 'stimnameL'] = stimnameL
        eventlist.at[eventnr, 'stimnameR'] = stimnameR
        eventlist.at[eventnr, 'stimnameC'] = stimnameC
        eventlist.at[eventnr, 'response'] = int(response)
        eventlist.at[eventnr, 'rt'] = rt-last_target if correct else 0
        eventlist.at[eventnr, 'correct'] = int(correct)
        eventlist.at[eventnr, 'time_stimon'] = time_stimon
        eventlist.at[eventnr, 'time_stimoff'] = time_stimoff
        eventlist.at[eventnr, 'stimdur'] = time_stimoff-time_stimon;
        
        while core.getTime() < time_stimon + isiduration:pass

finally:
    core.wait(1)
    fixation.draw()
    time_stimoff=win.flip()
    writeout(eventlist)
    progresstext.text = 'Experiment finished!'
    progresstext.draw()
    win.flip()
    core.wait(1)
    win.close()
    exit()


