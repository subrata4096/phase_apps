export LOOPITERSKIP_BOXBLUR='1'
export LOOPCHANGE_BOXBLUR='1'
export LOOPITERSKIP_NEIGHBOUR='1'
export LOOPCHANGE_NEIGHBOUR='1'
export LOOPCHANGE_HSCALE='1'
export LOOPITERSKIP_HSCALE='1'
export LOOPCHANGE_SWSCALE='1'
export LOOPITERSKIP_SWSCALE='1'
export LOOPCHANGE_H815SCALE='1'
export LOOPITERSKIP_H815SCALE='1'
export LOOPCHANGE_CONVOLVE='1'
export LOOPITERSKIP_CONVOLVE='1'
export TOTAL_FRAMES_MANISH='125'
export PHASETOAPPROXIMATE='-1'
export TOTAL_PHASES_MANISH='1'
./build/bin/ffmpeg -i videos/output/correct_out_def_blur_scale.mp4 -i videos/output/approxout_def_blur_scale.mp4 -lavfi  "ssim;[0:v][1:v]psnr" -f null -

./build/bin/ffmpeg -i ./videos/short_60fps_video_05ss.mp4 -vf deflate,boxblur=10,scale=1080:1920 ./videos/output/approxout_def_blur_scale.mp4 
./build/bin/ffmpeg -i ./videos/short_60fps_video_05ss.mp4 -vf deflate,boxblur=10,convolution="-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2:-2 -1 0 -1 1 1 0 1 2"  ./videos/output/approxout_def_blur_scale.mp4 
 ./build/bin/ffmpeg -i videos/top_doubt.mp4 -ss 00:00:00 -t 0:00:05 -async 1 ./videos/top_static_short.mp4

  ./build/bin/ffprobe -v error -count_frames -select_streams v:0   -show_entries stream=nb_read_frames -of default=nokey=1:noprint_wrappers=1 ./videos/adele_static_short.mp4 