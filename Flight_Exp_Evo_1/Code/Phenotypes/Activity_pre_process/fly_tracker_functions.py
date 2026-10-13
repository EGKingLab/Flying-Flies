import numpy as np
import cv2 as cv
import os
from pathlib import Path

def rmdir(directory):
    """
    Recursively remove a directory and its contents.

    Parameters:
    directory (str): Path to directory to remove.

    """
    directory = Path(directory)
    for item in directory.iterdir():
        if item.is_dir():
            rmdir(item)
        else:
            item.unlink()
    directory.rmdir()


def cleanup_directory(cleanup_dir):
    """
    Remove all splitvideos directories, labeled points movies, plate
    segmentation images, h5 files, and pickle files from a directory.

    Parameters:
    cleanup_dir (str): Path to directory to clean up.

    """

    for root, dirs, files in os.walk(cleanup_dir):
        # splitvideos folders
        print("Removing splitvideos directories")
        for name in dirs:
            if name.endswith('splitvideos'):
                print(os.path.join(root, name))
                rmdir(os.path.join(root, name))

        # labeled, segmentation, h5. pickle files
        print("Removing labeled points movies")
        for name in files:
            if name.endswith('_labeled.mp4'):
                print(os.path.join(root, name))
                os.remove(os.path.join(root, name))
        print("Removing plate segmentation images")
        for name in files:
            if name.endswith('_segmentation.jpg'):
                print(os.path.join(root, name))
                os.remove(os.path.join(root, name))
        print("Removing h5 files")
        for name in files:
            if name.endswith('.h5'):
                print(os.path.join(root, name))
                os.remove(os.path.join(root, name))
        print('\n\n')


def rotate_image(img, angle):
    """
    Rotate an image (np.array) by some arbitrary angle about the center of the image.

    See: https://stackoverflow.com/a/9042907/168137

    Parameters:
    img (np.array): Image as np.array, e.g., from opencv read()
    angle (float): Angle for rotation, negative angles are clockwise.

    Returns:
    np.array: Image of same dimensions as image.

    """

    image_center = tuple(np.array(img.shape[1::-1]) / 2)
    rot_mat = cv.getRotationMatrix2D(image_center, angle, 1.0)
    result = cv.warpAffine(img, rot_mat, img.shape[1::-1], flags=cv.INTER_LINEAR)
    return result


def crop(img, x1, y1, x2, y2):
    """
    Crop an image using the upper-left and lower-right corners of a rectangle.

    Parameters:
    img (np.array): Image as np.array, e.g., from opencv read()
    x1, y1 (int): Upper-left corner
    x2, y2 (int): Lower-right corner

    Returns:
    np.array: Cropped image.

    """

    img_crop = img[y1:y2, x1:x2, :]
    return img_crop


# https://www.pyimagesearch.com/2014/08/25/4-point-opencv-getperspective-transform-example/
def order_points(pts):
    """
    Order the points in a list of 4 points such that the first entry in the
    list is the top-left, the second entry is the top-right, the third is the
    bottom-right, and the fourth is the bottom-left.

    Parameters:
    pts (np.array): List of 4 points.

    """

    # initialzie a list of coordinates that will be ordered
    # such that the first entry in the list is the top-left,
    # the second entry is the top-right, the third is the
    # bottom-right, and the fourth is the bottom-left
    rect = np.zeros((4, 2), dtype = "float32")

    # the top-left point will have the smallest sum, whereas
    # the bottom-right point will have the largest sum
    s = pts.sum(axis = 1)
    rect[0] = pts[np.argmin(s)]
    rect[2] = pts[np.argmax(s)]

    # now, compute the difference between the points, the
    # top-right point will have the smallest difference,
    # whereas the bottom-left will have the largest difference
    diff = np.diff(pts, axis = 1)
    rect[1] = pts[np.argmin(diff)]
    rect[3] = pts[np.argmax(diff)]

    # return the ordered coordinates
    return rect


def four_point_transform(image, pts):
    """
    Perform a perspective transform on an image given 4 points.

    Parameters:
    image (np.array): Image as np.array, e.g., from opencv read()
    pts (np.array): List of 4 points.

    """

    # obtain a consistent order of the points and unpack them
    # individually
    rect = order_points(pts)
    (tl, tr, br, bl) = rect

    # compute the width of the new image, which will be the
    # maximum distance between bottom-right and bottom-left
    # x-coordiates or the top-right and top-left x-coordinates
    widthA = np.sqrt(((br[0] - bl[0]) ** 2) + ((br[1] - bl[1]) ** 2))
    widthB = np.sqrt(((tr[0] - tl[0]) ** 2) + ((tr[1] - tl[1]) ** 2))
    maxWidth = max(int(widthA), int(widthB))

    # compute the height of the new image, which will be the
    # maximum distance between the top-right and bottom-right
    # y-coordinates or the top-left and bottom-left y-coordinates
    heightA = np.sqrt(((tr[0] - br[0]) ** 2) + ((tr[1] - br[1]) ** 2))
    heightB = np.sqrt(((tl[0] - bl[0]) ** 2) + ((tl[1] - bl[1]) ** 2))
    maxHeight = max(int(heightA), int(heightB))

    # now that we have the dimensions of the new image, construct
    # the set of destination points to obtain a "birds eye view",
    # (i.e. top-down view) of the image, again specifying points
    # in the top-left, top-right, bottom-right, and bottom-left
    # order
    dst = np.array([
        [0, 0],
        [maxWidth - 1, 0],
        [maxWidth - 1, maxHeight - 1],
        [0, maxHeight - 1]], dtype = "float32")

    # compute the perspective transform matrix and then apply it
    M = cv.getPerspectiveTransform(rect, dst)
    warped = cv.warpPerspective(image, M, (maxWidth, maxHeight))

    # return the warped image
    return warped
