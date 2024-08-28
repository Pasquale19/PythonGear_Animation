import os
import re

# Directory containing the images
directory = "export"

# Function to extract number from filename
def get_number(filename):
    return int(re.search(r'\d+', filename).group())

# Get all image files and sort them numerically
image_files = [f for f in os.listdir(directory) if f.lower().endswith(('.png', '.jpg', '.jpeg', '.gif'))]
image_files.sort(key=get_number)

# Rename files
for i, old_name in enumerate(image_files, start=1):
    # Get file extension
    file_extension = os.path.splitext(old_name)[1]
    
    # Create new name with zero-padded number
    new_name = f"frame{i:03d}{file_extension}"
    
    # Full paths
    old_path = os.path.join(directory, old_name)
    new_path = os.path.join(directory, new_name)
    
    # Rename file
    os.rename(old_path, new_path)
    print(f"Renamed: {old_name} -> {new_name}")

print("Renaming complete!")