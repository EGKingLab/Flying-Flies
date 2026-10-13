#!/usr/bin/env python
# coding: utf-8

# ## Setup ####################################################################
import deeplabcut
import cv2 as cv
import os
import numpy as np
import time
import pandas as pd

from fly_tracker_functions import *

# ## Flags ####################################################################
write_tracked_video = False
cleanup_file_names = True

# ## Paths ####################################################################

# path to date directories to analyze
date_folders_base =  'D:\Fritz\DGRP\RawData'

###############################################################################
# Run this to check if folders are correct

# Option 1: Get list of all dates in date_folders_base (or supply specific 
# list below)
days_list = []
for file in os.listdir(date_folders_base):
    d = os.path.join(date_folders_base, file)
    if os.path.isdir(d):
        if file.startswith('20'):
            days_list.append(file)
    days_list.sort()

# Option 2: Provide a list of dates
# days_list = ['2021-12-08','2021-12-10', '2021-12-11', '2021-12-12',
#              '2021-12-13', '2021-12-14', '2021-12-15', '2021-12-16',
#              '2021-12-17', '2021-12-19', '2021-12-20', '2021-12-22']

print(days_list)

###############################################################################
###############################################################################
## Dont edit below this.
###############################################################################
###############################################################################

# Base path for config files
config_base_path = 'D:\FlyTracker'
os.chdir(config_base_path)

# deeplabcut fly tracker config.yaml path
config_path = os.path.join(config_base_path,
                           'fly_tracker_2-king-2021-09-27/config.yaml')

# deeplabcut plate corner finder config.yaml path
plate_config_path = os.path.join(config_base_path,
                                 'plate_finder-kinglab-2021-05-02/config.yaml')

# DLC model string
dlc_model_string = 'DLC_resnet50_fly_tracker_2Sep27shuffle1_1000000'

# Plate finder string
dlc_plate_string = 'DLC_resnet50_plate_finderMay2shuffle1_500000'

###############################################################################
###############################################################################

# Join the days_list with the date folder base to make full paths to the
# directories to be analyzed.
days_recording = [os.path.join(date_folders_base, x) for x in days_list]
print(days_recording)

# Plate coordinates
pos = list(range(0, 34))
x1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600,
      600, 600, 600]
y1 = [24, 79, 133, 187, 241, 297, 355, 413, 469, 524, 577, 634, 689, 744,
      799, 855, 908, 24, 79, 133, 187, 241, 297, 355, 413, 469, 524, 577,
      634, 689, 744, 799, 855, 908]
x2 = [599, 599, 599, 599, 599, 599, 599, 599, 599, 599, 599, 599, 599, 599,
      599, 599, 599, 1199, 1199, 1199, 1199, 1199, 1199, 1199, 1199, 1199,
      1199, 1199, 1199, 1199, 1199, 1199, 1199, 1199]
y2 = [78, 132, 186, 240, 296, 354, 412, 468, 523, 576, 633, 688, 743, 798,
      854, 907, 968, 78, 132, 186, 240, 296, 354, 412, 468, 523, 576, 633,
      688, 743, 798, 854, 907, 968]

pl = pd.DataFrame({'pos': pos, 'x1': x1, 'y1': y1, 'x2': x2, 'y2': y2})
pl['h'] = pl['y2'] - pl['y1']
pl['w'] = pl['x2'] - pl['x1']

# Iterate through folders
for jj, date in enumerate(days_recording):
    print("Processing " + date +"\n")

    # Cleanup directory
    print("Cleaning up directory")
    cleanup_directory(date)

    # Find .avi files
    files = []
    for r, d, f in os.walk(date):
        f.sort()
        for file in f:
            if file.endswith(".avi"):
                files.append(file)
    print(files)

    for ii, vidFile in enumerate(files):
        # Setup path for split video files
        splVidPath = os.path.join(date_folders_base,
                                  date,
                                  vidFile[:-4] + '_splitvideos')
        print('Videos will write into: ' + splVidPath)

        # Make directory for split videos
        try:
            os.mkdir(splVidPath)
        except OSError as error:
            print(error)

        # Find corners of plate
        deeplabcut.analyze_videos(plate_config_path,
                                  os.path.join(date_folders_base,
                                               date,
                                               vidFile),
                                  videotype='.avi',
                                  save_as_csv=False)

        # Write movie with plate corners
        time.sleep(10)
        deeplabcut.create_labeled_video(plate_config_path,
                                  os.path.join(date_folders_base,
                                               date,
                                               vidFile),
                                  videotype='.avi')

        # Rename plate corner finder output files
        old_file =  os.path.join(date_folders_base,
                                 date, vidFile[:-4] + dlc_plate_string + '.h5')
        os.rename(old_file, old_file.replace(dlc_plate_string, ''))

        old_file = os.path.join(date_folders_base,
                                date, vidFile[:-4] + dlc_plate_string + '_labeled.mp4')
        os.rename(old_file, old_file.replace(dlc_plate_string, ''))

        # Remove pickle created by DLC
        os.remove(os.path.join(date_folders_base,
                               date, vidFile[:-4] + dlc_plate_string + '_meta.pickle'))

        # Create videoWriters
        for vw_id in range(34):
            outfile = os.path.join(splVidPath,
                                   vidFile[:-4] + '_' + str(vw_id) + '.mp4')
            globals()['vW{}'.format(vw_id)] = cv.VideoWriter(
                filename=outfile,
                apiPreference=cv.CAP_FFMPEG,
                fourcc=cv.VideoWriter_fourcc('m', 'p', '4', 'v'),
                fps=4,
                frameSize=(pl.w[vw_id], pl.h[vw_id]))


        # Processing video ###############################################
        # Load video
        vid = cv.VideoCapture(os.path.join(date_folders_base, date, vidFile))
        nframes = int(vid.get(cv.CAP_PROP_FRAME_COUNT))

        # Load csv with tracked corners
        corners = pd.read_hdf(os.path.join(date_folders_base, date,
                              vidFile[:-4] + '.h5'))
        corners = corners[dlc_plate_string]

        # Extract coordinates
        UL = corners['topleft']
        UR = corners['topright']
        LL = corners['bottomleft']
        LR = corners['bottomright']

        # Iterate through frames
        print("\nSplitting into separate videos and cropping.\n")
        for ii in range(nframes):

            # Load image
            ret, img = vid.read()

            # Corner points
            corner_pts = np.array([
                [LL.iloc[ii, 0], LL.iloc[ii, 1]],
                [UL.iloc[ii, 0], UL.iloc[ii, 1]],
                [UR.iloc[ii, 0], UR.iloc[ii, 1]],
                [LR.iloc[ii, 0], LR.iloc[ii, 1]]
                ])

            # Use four point transformation
            warped = four_point_transform(img, corner_pts)

            # Resize to 1200 x 1000
            dim = (1200, 1000)

            # resize image
            resized = cv.resize(warped, dim, interpolation = cv.INTER_AREA)

            # Plate segmentation check. Use 100th frame.
            if ii == 99:
                seg_check = resized.copy()
                for idx, row in pl.iterrows():
                    cv.rectangle(seg_check,
                                 (row.x1, row.y1),
                                 (row.x2, row.y2),
                                 (255,0,0), 1)
                    cv.putText(seg_check, str(idx), (row.x1+15, row.y1+15),
                               cv.FONT_HERSHEY_SIMPLEX, 0.5, (255,0,0), 1, cv.LINE_AA)
                ret = cv.imwrite(os.path.join(date_folders_base,
                                              date,
                                              vidFile[:-4] + '_segmentation.jpg'), seg_check)

            ## Segment video and write frame to that videoWriter
            for idx, row in pl.iterrows():
                # Create copy of img
                img_tmp = resized.copy()

                # Crop image
                img_crop = crop(img_tmp, row.x1, row.y1, row.x2, row.y2)
                h,w = img_crop.shape[:2]

                # Write image
                globals()['vW{}'.format(idx)].write(img_crop)

        # Release videoWriters
        for vw_id in range(34):
            globals()['vW{}'.format(vw_id)].release()

        # Release main video file
        vid.release()

        # Track points on cropped videos
        print('Pausing before tracking points.\n')
        time.sleep(10)
        deeplabcut.analyze_videos(config_path,
                                  [splVidPath],
                                  videotype='.mp4',
                                  save_as_csv=True)

        # Write movie with tracked fly
        if write_tracked_video:
            print('\nPausing before writing videos with tracked points.\n')
            time.sleep(10)
            deeplabcut.create_labeled_video(config_path,
                                            [splVidPath],
                                            videotype='.mp4',
                                            trailpoints=20)

        # rename output files
        if cleanup_file_names:
            for r, d, f in os.walk(splVidPath):
                for file in f:
                    if file != None and dlc_model_string in file:
                        os.rename(os.path.join(r, file),
                                  os.path.join(r, file.replace(dlc_model_string, '')))

        # Pause before starting the next video
        print('\nPausing before processing the next video.\n')
        time.sleep(20)
