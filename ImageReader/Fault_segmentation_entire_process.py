import cv2
import numpy as np
import math
from skimage.morphology import skeletonize
from skimage.util import invert
import sknw
import networkx as nx
from scipy.ndimage import convolve
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra


"""
Input file name should be [Faultimage.jpg]. 

This script processes an input drawing of a fault line by performing the following steps:
1. Background removal: Keeps only the black fault lines and green scale line parts of the image.
2. Fault segmentation: Removes the green scale line, skeletonizes the remaining drawing,
   computes its length using Dijkstra's algorithm, and segments it into individual lines.
3. Line segment splitting: Filters out duplicate or very similar overlapping segments.

The processed outputs include images with marked endpoints and a text file with the segment endpoints,
with coordinates scaled to real-world units (kilometers).
"""


# Give the scale length (units correspond to its real-life length in kilometers)
scale_length = 10.0

# Define the min and max segment lengths in real-world units, then convert to pixels
min_length_units = 2.0  # minimum segment length in units
max_length_units = 5.0   # maximum segment length in units

# Change this line to rename the file you want to analyze
base_name = "ImageReader/Faultimage"

# Write the real file names and paths
drawing_file = f"{base_name}.jpg"
no_background_file = f"{base_name}_drawing_no_background.jpg"
segmented_endpoints_file = f"{base_name}_segmented_endpoints.txt"
segmented_file = f"{base_name}_segmented.jpg"


def background_remover():
    """
    Reads the drawing file, removes the background by keeping only black and green parts,
    and saves the result as the no-background file.
    """
    # Read the image
    image = cv2.imread(drawing_file)
    if image is None:
        print(f"Can't read the image '{drawing_file}'.")
        exit()

    # Convert to HSV color space
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)

    # Define HSV ranges for black and green
    lower_black = np.array([0, 0, 0])
    upper_black = np.array([180, 255, 50])
    lower_green = np.array([35, 40, 40])
    upper_green = np.array([85, 255, 255])

    # Create masks for black and green colors
    mask_black = cv2.inRange(hsv, lower_black, upper_black)
    mask_green = cv2.inRange(hsv, lower_green, upper_green)

    # Combine the masks
    mask = cv2.bitwise_or(mask_black, mask_green)

    # Create an output image with a white background
    output = np.full_like(image, 255)
    output[mask != 0] = image[mask != 0]

    # Save the image without background
    cv2.imwrite(no_background_file, output)
    print("Background removal done")

def fault_segmentation():
    """
    Reads the no-background image, takes out the green fault line,
    computes a skeleton and its length, then segments the skeleton 
    into lines, and finally saves both an image with blue points
    and a text file with the endpoints.
    """
    # Read the image with no background
    image = cv2.imread(no_background_file)
    if image is None:
        print(f"Can't read '{no_background_file}'.")
        exit()
    height, width = image.shape[:2]

    # Convert image to HSV and create a mask for green
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)
    lower_green = np.array([40, 40, 40])
    upper_green = np.array([80, 255, 255])
    mask_green = cv2.inRange(hsv, lower_green, upper_green)
    mask_green_bool = mask_green > 0

    # Skeletonize the green mask
    skeleton_green = skeletonize(mask_green_bool)

    def find_endpoints(skel):
        """
        Finds endpoints in a skeletonized binary image.
        """
        kernel = np.array([[1, 1, 1],
                           [1, 10, 1],
                           [1, 1, 1]])
        neighbor_count = convolve(skel.astype(int), kernel, mode='constant', cval=0)
        endpoints = np.logical_and(skel, neighbor_count == 11)
        return endpoints

    # Find endpoints of the green skeleton and mark them with blue dots
    endpoints_green = find_endpoints(skeleton_green)
    coords = np.column_stack(np.where(endpoints_green))
    if coords.shape[0] != 2:
        print("Found more than two endpoints for the green line.")
    # In OpenCV, color is in BGR. Blue is (255, 0, 0)
    for y, x in coords:
        cv2.circle(image, (x, y), radius=3, color=(255, 0, 0), thickness=-1)

    def compute_skeleton_length(skel):
        """
        Calculates the length of a skeleton (in pixels) using a Dijkstra shortest-path
        between its two endpoints.
        """
        coords = np.column_stack(np.where(skel))
        indices = np.ravel_multi_index((coords[:, 0], coords[:, 1]), skel.shape)
        index_map = -np.ones(skel.shape, dtype=int)
        index_map[coords[:, 0], coords[:, 1]] = np.arange(len(coords))
        data = []
        row_ind = []
        col_ind = []
        for i in range(len(coords)):
            y, x = coords[i]
            neighbors = [(y - 1, x - 1), (y - 1, x), (y - 1, x + 1),
                         (y, x - 1),               (y, x + 1),
                         (y + 1, x - 1), (y + 1, x), (y + 1, x + 1)]
            for ny, nx in neighbors:
                if 0 <= ny < skel.shape[0] and 0 <= nx < skel.shape[1]:
                    if skel[ny, nx]:
                        j = index_map[ny, nx]
                        if j >= 0:
                            dist = np.hypot(ny - y, nx - x)
                            data.append(dist)
                            row_ind.append(i)
                            col_ind.append(j)
        adjacency = csr_matrix((data, (row_ind, col_ind)), shape=(len(coords), len(coords)))
        endpoints_coords = np.column_stack(np.where(find_endpoints(skel)))
        if endpoints_coords.shape[0] != 2:
            print("Warning: Skeleton doesn't have exactly two endpoints.")
            return None
        start_idx = index_map[endpoints_coords[0, 0], endpoints_coords[0, 1]]
        end_idx = index_map[endpoints_coords[1, 0], endpoints_coords[1, 1]]
        dist_matrix, predecessors = dijkstra(csgraph=adjacency, directed=False,
                                             indices=start_idx, return_predecessors=True)
        length = dist_matrix[end_idx]
        return length

    # Compute the length of the green line and the scale factor (units correspond to its real-life length in kilometers)
    length_in_pixels = compute_skeleton_length(skeleton_green)
    if length_in_pixels is None:
        print("Can't compute the length of the green line.")
        exit()
    scale_factor = scale_length / length_in_pixels

    # Remove the green parts (fault line) from the image by setting those pixels to white
    image[mask_green_bool] = [255, 255, 255]

    for y, x in coords:
        cv2.circle(image, (x, y), radius=4, color=(255, 255, 255), thickness=-1)

    # Prepare a binary image for skeletonization of the remaining (fault) drawing
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    _, binary = cv2.threshold(gray, 50, 255, cv2.THRESH_BINARY_INV)
    binary_bool = binary > 0
    skeleton = skeletonize(binary_bool)

    # Build a graph from the skeleton using sknw
    graph = sknw.build_sknw(skeleton.astype(np.uint8))

    min_length_pixels = min_length_units / scale_factor
    max_length_pixels = max_length_units / scale_factor

    def segment_edge(pts, min_len, max_len):
        """
        Splits a sequence of points (edge) into segments whose length is between min_len and max_len.
        """
        dists = np.sqrt(np.sum(np.diff(pts, axis=0) ** 2, axis=1))
        cumulative_dists = np.concatenate(([0], np.cumsum(dists)))
        segments = []
        start_idx = 0
        while start_idx < len(pts) - 1:
            end_idx = start_idx + 1
            while end_idx < len(pts):
                delta = cumulative_dists[end_idx] - cumulative_dists[start_idx]
                if delta >= min_len and delta <= max_len:
                    segments.append((pts[start_idx], pts[end_idx]))
                    start_idx = end_idx
                    break
                elif delta > max_len:
                    segments.append((pts[start_idx], pts[end_idx - 1]))
                    start_idx = end_idx - 1
                    break
                else:
                    end_idx += 1
            else:
                segments.append((pts[start_idx], pts[-1]))
                break
        return segments

    # Process each edge in the graph to split it into segments
    line_segments = []
    for (s, e), edge in graph.edges.items():
        pts = edge['pts']
        segments = segment_edge(pts, min_length_pixels, max_length_pixels)
        line_segments.extend(segments)

    # Save the endpoints of each segment into a text file.
    # The coordinates are transformed to a coordinate system where the origin is at the bottom left,
    # and scaled by the previously found scale factor.
    with open(segmented_endpoints_file, 'w') as f:
        for seg in line_segments:
            p1, p2 = seg
            y1_img, x1_img = p1
            y2_img, x2_img = p2

            # Convert to a standard coordinate system (with the origin at bottom left)
            y1_std = height - y1_img
            y2_std = height - y2_img

            # Apply the scale factor
            x1 = x1_img * scale_factor
            x2 = x2_img * scale_factor
            y1 = y1_std * scale_factor
            y2 = y2_std * scale_factor

            f.write(f"{x1}, {y1}, {x2}, {y2}\n")

    # Mark the endpoints of each segment on the image with blue dots (BGR: (255, 0, 0))
    for seg in line_segments:
        p1, p2 = seg
        y1, x1 = p1
        y2, x2 = p2
        cv2.circle(image, (int(x1), int(y1)), radius=2, color=(255, 0, 0), thickness=-1)
        cv2.circle(image, (int(x2), int(y2)), radius=2, color=(255, 0, 0), thickness=-1)

    # Save the segmented image with marked endpoints
    cv2.imwrite(segmented_file, image)
    print("Fault segmentation done")

def line_segment_splitter():
    """
    Reads the endpoints file, removes segments that are too close to each other,
    and writes the filtered endpoints back to the same file.
    """
    rows = []
    with open(segmented_endpoints_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            # Split by comma and convert the strings to floats
            a_str, b_str, c_str, d_str = line.split(',')
            a = float(a_str)
            b = float(b_str)
            c = float(c_str)
            d = float(d_str)
            rows.append((a, b, c, d))

    # Identify rows that are very similar (like duplicate segments) and delete them
    rows_to_delete = set()
    for i in range(len(rows)):
        if i in rows_to_delete:
            continue
        a_i, b_i, c_i, d_i = rows[i]
        for j in range(i + 1, len(rows)):
            if j in rows_to_delete:
                continue
            a_j, b_j, c_j, d_j = rows[j]
            dist1 = math.hypot(a_i - a_j, b_i - b_j)
            dist2 = math.hypot(c_i - c_j, d_i - d_j)
            if dist1 < 0.01 and dist2 < 0.01:
                rows_to_delete.add(i)
                break

    # Write back only the unique rows
    with open(segmented_endpoints_file, 'w') as f_new:
        for idx, row in enumerate(rows):
            if idx not in rows_to_delete:
                f_new.write(', '.join(map(str, row)) + '\n')

    print("Line segment splitting done")

if __name__ == "__main__":
    print("Background removal")
    background_remover()

    print("Fault segmentation")
    fault_segmentation()

    print("Line segment split")
    line_segment_splitter()

    print("All steps done")