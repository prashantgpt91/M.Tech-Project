import numpy as np
import cv2

MIN_MATCH_COUNT = 50
object_img = cv2.imread('13.jpg',0) # queryImage
scene_img = cv2.imread('12.jpg',0) # trainImage

def checkduplicate(object_img, scene_img):
    # Initiate SIFT detector
    akaze = cv2.AKAZE_create()

    # find the keypoints and descriptors with SIFT
    kp1, des1 = akaze.detectAndCompute(object_img,None)
    kp2, des2 = akaze.detectAndCompute(scene_img,None)

    FLANN_INDEX_KDTREE = 0

    bf = cv2.BFMatcher()
    matches = bf.knnMatch(des1,des2, k=2)

    # store all the good matches as per Lowe's ratio test.
    good = []
    for m,n in matches:
        if m.distance < 0.7*n.distance:
            good.append(m)
    print len(good)
    if len(good)>MIN_MATCH_COUNT:
        src_pts = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1,1,2)
        dst_pts = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1,1,2)

        M, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC,5.0)
        h,w = object_img.shape
        pts = np.float32([ [0,0],[0,h-1],[w-1,h-1],[w-1,0] ]).reshape(-1,1,2)
        dst = cv2.perspectiveTransform(pts,M)
        a = dst.astype(int)
        a = a.flatten()
        a = a.reshape(4, 2)
        crop = scene_img[a[0][0]:a[0][0] + h, a[0][1]:a[0][1] + w]
        cv2.imshow("outjpg",crop)
        return crop

    else:
        print "Not enough matches are found - %d/%d" % (len(good),MIN_MATCH_COUNT)
        return None

img = checkduplicate(object_img,scene_img)
cv2.waitKey(0)