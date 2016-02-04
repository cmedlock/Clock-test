#!/bin/bash

DATADIR=~/Documents/DSP_UROP/all_data

# create directories
mkdir -p $DATADIR/animations/command_healthy
mkdir -p $DATADIR/animations/command_impaired
mkdir -p $DATADIR/animations/copy_healthy
mkdir -p $DATADIR/animations/copy_impaired

mkdir -p $DATADIR/animations_3d/command_healthy
mkdir -p $DATADIR/animations_3d/command_impaired
mkdir -p $DATADIR/animations_3d/copy_healthy
mkdir -p $DATADIR/animations_3d/copy_impaired

mkdir -p $DATADIR/xy_true/command_healthy
mkdir -p $DATADIR/xy_true/command_impaired
mkdir -p $DATADIR/xy_true/copy_healthy
mkdir -p $DATADIR/xy_true/copy_impaired

mkdir -p $DATADIR/xy_polar/command_healthy
mkdir -p $DATADIR/xy_polar/command_impaired
mkdir -p $DATADIR/xy_polar/copy_healthy
mkdir -p $DATADIR/xy_polar/copy_impaired

mkdir -p $DATADIR/velocity_acceleration_xy/command_healthy
mkdir -p $DATADIR/velocity_acceleration_xy/command_impaired
mkdir -p $DATADIR/velocity_acceleration_xy/copy_healthy
mkdir -p $DATADIR/velocity_acceleration_xy/copy_impaired

mkdir -p $DATADIR/velocity_acceleration_overall/command_healthy
mkdir -p $DATADIR/velocity_acceleration_overall/command_impaired
mkdir -p $DATADIR/velocity_acceleration_overall/copy_healthy
mkdir -p $DATADIR/velocity_acceleration_overall/copy_impaired

mkdir -p $DATADIR/pressure/command_healthy
mkdir -p $DATADIR/pressure/command_impaired
mkdir -p $DATADIR/pressure/copy_healthy
mkdir -p $DATADIR/pressure/copy_impaired

# transfer files
#cp $DATADIR/figs_raw/YDU*/command_clock_anim*mp4 $DATADIR/animations/command_healthy
#cp $DATADIR/figs_raw/CIN*/command_clock_anim*mp4 $DATADIR/animations/command_impaired
#cp $DATADIR/figs_raw/YDU*/copy_clock_anim*mp4 $DATADIR/animations/copy_healthy
#cp $DATADIR/figs_raw/CIN*/copy_clock_anim*mp4 $DATADIR/animations/copy_impaired

#cp $DATADIR/figs_raw/YDU*/COMMAND_clock_3danim*mp4 $DATADIR/animations_3d/command_healthy
#cp $DATADIR/figs_raw/CIN*/COMMAND_clock_3danim*mp4 $DATADIR/animations_3d/command_impaired
#cp $DATADIR/figs_raw/YDU*/COPY_clock_3danim*mp4 $DATADIR/animations_3d/copy_healthy
#cp $DATADIR/figs_raw/CIN*/COPY_clock_3danim*mp4 $DATADIR/animations_3d/copy_impaired

#cp $DATADIR/figs_raw/YDU*/whole_COMMAND*png $DATADIR/animations_3d/command_healthy
#cp $DATADIR/figs_raw/CIN*/whole_COMMAND*png $DATADIR/animations_3d/command_impaired
#cp $DATADIR/figs_raw/YDU*/whole_COPY*png $DATADIR/animations_3d/copy_healthy
#cp $DATADIR/figs_raw/CIN*/whole_COPY*png $DATADIR/animations_3d/copy_impaired

#cp $DATADIR/figs_raw/YDU*/*_true_command*png $DATADIR/xy_true/command_healthy
#cp $DATADIR/figs_raw/CIN*/*_true_command*png $DATADIR/xy_true/command_impaired
#cp $DATADIR/figs_raw/YDU*/*_true_copy*png $DATADIR/xy_true/copy_healthy
#cp $DATADIR/figs_raw/CIN*/*_true_copy*png $DATADIR/xy_true/copy_impaired

#cp $DATADIR/figs_raw/YDU*/xy_command*png $DATADIR/xy_true/command_healthy
#cp $DATADIR/figs_raw/CIN*/xy_command*png $DATADIR/xy_true/command_impaired
#cp $DATADIR/figs_raw/YDU*/xy_copy*png $DATADIR/xy_true/copy_healthy
#cp $DATADIR/figs_raw/CIN*/xy_copy*png $DATADIR/xy_true/copy_impaired

cp $DATADIR/figs_raw/YDU*/xy_polar_COMMAND*png $DATADIR/xy_polar/command_healthy
cp $DATADIR/figs_raw/CIN*/xy_polar_COMMAND*png $DATADIR/xy_polar/command_impaired
cp $DATADIR/figs_raw/YDU*/xy_polar_COPY*png $DATADIR/xy_polar/copy_healthy
cp $DATADIR/figs_raw/CIN*/xy_polar_COPY*png $DATADIR/xy_polar/copy_impaired

#cp $DATADIR/figs_raw/YDU*/v*t_a*t_command*png $DATADIR/velocity_acceleration_xy/command_healthy
#cp $DATADIR/figs_raw/CIN*/v*t_a*t_command*png $DATADIR/velocity_acceleration_xy/command_impaired
#cp $DATADIR/figs_raw/YDU*/v*t_a*t_copy*png $DATADIR/velocity_acceleration_xy/copy_healthy
#cp $DATADIR/figs_raw/CIN*/v*t_a*t_copy*png $DATADIR/velocity_acceleration_xy/copy_impaired

#cp $DATADIR/figs_raw/YDU*/vt_at_command*png $DATADIR/velocity_acceleration_overall/command_healthy
#cp $DATADIR/figs_raw/CIN*/vt_at_command*png $DATADIR/velocity_acceleration_overall/command_impaired
#cp $DATADIR/figs_raw/YDU*/vt_at_copy*png $DATADIR/velocity_acceleration_overall/copy_healthy
#cp $DATADIR/figs_raw/CIN*/vt_at_copy*png $DATADIR/velocity_acceleration_overall/copy_impaired

#cp $DATADIR/figs_raw/YDU*/pt_command*png $DATADIR/pressure/command_healthy
#cp $DATADIR/figs_raw/CIN*/pt_command*png $DATADIR/pressure/command_impaired
#cp $DATADIR/figs_raw/YDU*/pt_copy*png $DATADIR/pressure/copy_healthy
#cp $DATADIR/figs_raw/CIN*/pt_copy*png $DATADIR/pressure/copy_impaired

echo Done
