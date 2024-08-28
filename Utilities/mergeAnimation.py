import os
import subprocess

#overwriteing not possible

# Define the list of MP4 files to merge
video_files = [
    "animation1.mp4",
    "animation2.mp4"
]

video_files = [
    "animation1-converted.mp4",
    "animation2-converted.mp4"
]

video_files = [    "animation1.mp4",    "animation2.mp4","animation3.mp4","animation4.mp4"]

# Define the output directory and output video file
output_directory = os.path.expanduser("~/Videos")
output_video = os.path.join(output_directory, "merged2.mp4")

# Create a text file with the list of video files
with open("mylist.txt", "w") as f:
    for video_file in video_files:
        f.write(f"file '{os.path.abspath(video_file)}'\n")

# Merge the videos using FFmpeg
ffmpeg_command = [
    "ffmpeg",
    "-f", "concat",
    "-safe", "0",
    "-i", "mylist.txt",
    "-c", "copy",
    "-v", "verbose",  # Add verbose output
    output_video
]

try:
    result = subprocess.run(ffmpeg_command, check=True, capture_output=True, text=True)
    print("FFmpeg Output:", result.stdout)
except subprocess.CalledProcessError as e:
    print("Error occurred:")
    print("FFmpeg Output:", e.stdout)
    print("FFmpeg Error:", e.stderr)

# Clean up the mylist.txt file
os.remove("mylist.txt")

if os.path.exists(output_video):
    print(f"Merged video saved as: {output_video}\n merged: {video_files}")
else:
    print("Failed to create merged video.")
print(f"finished overwrite not possible")
    #vlc command disable repeat
#vlc_location ["animation1.mp4"] ["animation2.mp4"] --sout "#gather:std{access=file,dst=[combined.mp4]}" --sout-keep
#"C:\Program Files\VideoLAN\VLC\vlc.exe" "animation1-converted.mp4" "animation2-converted.mp4" --sout "#gather:std{access=file,dst=combined.mp4}" --sout-keep

#"C:\Program Files\VideoLAN\VLC\vlc.exe" -I dummy "animation1-converted.mp4" "animation2-converted.mp4" --sout "#gather:std{access=file,mux=mp4,dst=combined2.mp4}" --no-sout-all --sout-keep vlc://quit

#ffprobe -v error -select_streams v:0 -show_entries stream=codec_name -of default=noprint_wrappers=1:nokey=1 animation1-converted.mp4