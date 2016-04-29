#!/bin/bash

DATADIR=~/Documents/DSP_UROP/DataCatherine
DROPBOX=~/Dropbox\ \(MIT\)/DSP\ UROP/

# create directories
#mkdir -p $DATADIR/clock_face_animations_2d/command_healthy
#mkdir -p $DATADIR/clock_face_animations_2d/command_impaired
#mkdir -p $DATADIR/clock_face_animations_2d/copy_healthy
#mkdir -p $DATADIR/clock_face_animations_2d/copy_impaired

#mkdir -p $DATADIR/DFS_coefficients/command_healthy
#mkdir -p $DATADIR/DFS_coefficients/command_impaired
#mkdir -p $DATADIR/DFS_coefficients/copy_healthy
#mkdir -p $DATADIR/DFS_coefficients/copy_impaired

#mkdir -p $DATADIR/whole_clock_animations_3d/command_healthy
#mkdir -p $DATADIR/whole_clock_animations_3d/command_impaired
#mkdir -p $DATADIR/whole_clock_animations_3d/copy_healthy
#mkdir -p $DATADIR/whole_clock_animations_3d/copy_impaired

#mkdir -p $DATADIR/clock_face_figures_2d/command_healthy
#mkdir -p $DATADIR/clock_face_figures_2d/command_impaired
#mkdir -p $DATADIR/clock_face_figures_2d/copy_healthy
#mkdir -p $DATADIR/clock_face_figures_2d/copy_impaired

#mkdir -p $DATADIR/whole_clock_figures_2d/command_healthy
#mkdir -p $DATADIR/whole_clock_figures_2d/command_impaired
#mkdir -p $DATADIR/whole_clock_figures_2d/copy_healthy
#mkdir -p $DATADIR/whole_clock_figures_2d/copy_impaired

#mkdir -p $DATADIR/whole_clock_figures_3d/command_healthy
#mkdir -p $DATADIR/whole_clock_figures_3d/command_impaired
#mkdir -p $DATADIR/whole_clock_figures_3d/copy_healthy
#mkdir -p $DATADIR/whole_clock_figures_3d/copy_impaired

#mkdir -p $DATADIR/xy_true/command_healthy
#mkdir -p $DATADIR/xy_true/command_impaired
#mkdir -p $DATADIR/xy_true/copy_healthy
#mkdir -p $DATADIR/xy_true/copy_impaired

#mkdir -p $DATADIR/xy_polar/command_healthy
#mkdir -p $DATADIR/xy_polar/command_impaired
#mkdir -p $DATADIR/xy_polar/copy_healthy
#mkdir -p $DATADIR/xy_polar/copy_impaired

#mkdir -p $DATADIR/LPC/command_healthy
#mkdir -p $DATADIR/LPC/command_impaired
#mkdir -p $DATADIR/LPC/copy_healthy
#mkdir -p $DATADIR/LPC/copy_impaired

mkdir -p $DATADIR/LPC_complex/command_healthy
mkdir -p $DATADIR/LPC_complex/command_impaired
mkdir -p $DATADIR/LPC_complex/copy_healthy
mkdir -p $DATADIR/LPC_complex/copy_impaired

#mkdir -p $DATADIR/velocity_acceleration_xy/command_healthy
#mkdir -p $DATADIR/velocity_acceleration_xy/command_impaired
#mkdir -p $DATADIR/velocity_acceleration_xy/copy_healthy
#mkdir -p $DATADIR/velocity_acceleration_xy/copy_impaired

#mkdir -p $DATADIR/velocity_acceleration_overall/command_healthy
#mkdir -p $DATADIR/velocity_acceleration_overall/command_impaired
#mkdir -p $DATADIR/velocity_acceleration_overall/copy_healthy
#mkdir -p $DATADIR/velocity_acceleration_overall/copy_impaired

#mkdir -p $DATADIR/pressure/command_healthy
#mkdir -p $DATADIR/pressure/command_impaired
#mkdir -p $DATADIR/pressure/copy_healthy
#mkdir -p $DATADIR/pressure/copy_impaired

# transfer files
#cp $DATADIR/figs_raw/YDU*/command_clock_anim*mp4 $DATADIR/clock_face_animations_2d/command_healthy
#cp $DATADIR/figs_raw/CIN*/command_clock_anim*mp4 $DATADIR/clock_face_animations_2d/command_impaired
#cp $DATADIR/figs_raw/YDU*/copy_clock_anim*mp4 $DATADIR/clock_face_animations_2d/copy_healthy
#cp $DATADIR/figs_raw/CIN*/copy_clock_anim*mp4 $DATADIR/clock_face_animations_2d/copy_impaired

#cp $DATADIR/figs_raw/YDU*/dft*_COMMAND*png $DATADIR/DFS_coefficients/command_healthy
#cp $DATADIR/figs_raw/CIN*/dft*_COMMAND*png $DATADIR/DFS_coefficients/command_impaired
#cp $DATADIR/figs_raw/YDU*/dft*_COPY*png $DATADIR/DFS_coefficients/copy_healthy
#cp $DATADIR/figs_raw/CIN*/dft*_COPY*png $DATADIR/DFS_coefficients/copy_impaired

#cp $DATADIR/figs_raw/YDU*/COMMAND_clock_3danim*mp4 $DATADIR/whole_clock_animations_3d/command_healthy
#cp $DATADIR/figs_raw/CIN*/COMMAND_clock_3danim*mp4 $DATADIR/whole_clock_animations_3d/command_impaired
#cp $DATADIR/figs_raw/YDU*/COPY_clock_3danim*mp4 $DATADIR/whole_clock_animations_3d/copy_healthy
#cp $DATADIR/figs_raw/CIN*/COPY_clock_3danim*mp4 $DATADIR/whole_clock_animations_3d/copy_impaired

#cp $DATADIR/figs_raw/YDU*/xy_COMMAND*png $DATADIR/clock_face_figures_2d/command_healthy
#cp $DATADIR/figs_raw/CIN*/xy_COMMAND*png $DATADIR/clock_face_figures_2d/command_impaired
#cp $DATADIR/figs_raw/YDU*/xy_COPY*png $DATADIR/clock_face_figures_2d/copy_healthy
#cp $DATADIR/figs_raw/CIN*/xy_COPY*png $DATADIR/clock_face_figures_2d/copy_impaired

#cp $DATADIR/figs_raw/YDU*/whole_COMMAND_clock_3d*png $DATADIR/whole_clock_figures_3d/command_healthy
#cp $DATADIR/figs_raw/CIN*/whole_COMMAND_clock_3d*png $DATADIR/whole_clock_figures_3d/command_impaired
#cp $DATADIR/figs_raw/YDU*/whole_COPY_clock_3d*png $DATADIR/whole_clock_figures_3d/copy_healthy
#cp $DATADIR/figs_raw/CIN*/whole_COPY_clock_3d*png $DATADIR/whole_clock_figures_3d/copy_impaired

#cp $DATADIR/figs_raw/YDU*/whole_COMMAND_clock_2d*png $DATADIR/whole_clock_figures_2d/command_healthy
#cp $DATADIR/figs_raw/CIN*/whole_COMMAND_clock_2d*png $DATADIR/whole_clock_figures_2d/command_impaired
#cp $DATADIR/figs_raw/YDU*/whole_COPY_clock_2d*png $DATADIR/whole_clock_figures_2d/copy_healthy
#cp $DATADIR/figs_raw/CIN*/whole_COPY_clock_2d*png $DATADIR/whole_clock_figures_2d/copy_impaired

#cp $DATADIR/figs_raw/YDU*/*_true_COMMAND*png $DATADIR/xy_true/command_healthy
#cp $DATADIR/figs_raw/CIN*/*_true_COMMAND*png $DATADIR/xy_true/command_impaired
#cp $DATADIR/figs_raw/YDU*/*_true_COPY*png $DATADIR/xy_true/copy_healthy
#cp $DATADIR/figs_raw/CIN*/*_true_COPY*png $DATADIR/xy_true/copy_impaired

#cp $DATADIR/figs_raw/YDU*/xy_polar_COMMAND*png $DATADIR/xy_polar/command_healthy
#cp $DATADIR/figs_raw/CIN*/xy_polar_COMMAND*png $DATADIR/xy_polar/command_impaired
#cp $DATADIR/figs_raw/YDU*/xy_polar_COPY*png $DATADIR/xy_polar/copy_healthy
#cp $DATADIR/figs_raw/CIN*/xy_polar_COPY*png $DATADIR/xy_polar/copy_impaired

#cp $DATADIR/figs_raw/YDU*/lpc*covFalse*COMMAND*png $DATADIR/LPC/command_healthy
#cp $DATADIR/figs_raw/CIN*/lpc*covFalse*COMMAND*png $DATADIR/LPC/command_impaired
#cp $DATADIR/figs_raw/YDU*/lpc*covFalse*COPY*png $DATADIR/LPC/copy_healthy
#cp $DATADIR/figs_raw/CIN*/lpc*covFalse*COPY*png $DATADIR/LPC/copy_impaired

cp $DATADIR/figs_raw/YDU*/complex_lpc*COMMAND*png $DATADIR/LPC_complex/command_healthy
cp $DATADIR/figs_raw/CIN*/complex_lpc*COMMAND*png $DATADIR/LPC_complex/command_impaired
cp $DATADIR/figs_raw/YDU*/complex_lpc*COPY*png $DATADIR/LPC_complex/copy_healthy
cp $DATADIR/figs_raw/CIN*/complex_lpc*COPY*png $DATADIR/LPC_complex/copy_impaired

#cp $DATADIR/figs_raw/YDU*/instfreq*COMMAND*png $DATADIR/LPC_complex/command_healthy
#cp $DATADIR/figs_raw/CIN*/instfreq*COMMAND*png $DATADIR/LPC_complex/command_impaired
#cp $DATADIR/figs_raw/YDU*/instfreq*COPY*png $DATADIR/LPC_complex/copy_healthy
#cp $DATADIR/figs_raw/CIN*/instfreq*COPY*png $DATADIR/LPC_complex/copy_impaired

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
