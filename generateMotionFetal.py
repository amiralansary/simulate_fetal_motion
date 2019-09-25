# library path
import os, sys
lib_path = os.path.abspath('/vol/biomedic/users/aa16914/lib/python2.7/site-packages/')
sys.path.insert(1,lib_path)
lib_path = os.path.abspath('/vol/biomedic/users/aa16914/lib/SimpleITK/build/Wrapping/build/lib.linux-x86_64-2.7/SimpleITK/')
sys.path.insert(1,lib_path)


import  SimpleITK   as sitk
import  numpy       as np
import  math        as mt


def resample_sitk(img, newSpacing, shiftOrigin=(0,0,0), interpolator=sitk.sitkBSpline):
    """ This function transforms a nifti image
        Attributes:
            img:            The fixed image that will be transformed (simpleitk type)
            newSpacing:     The translation vector in mm [tx,ty,tz]
            shiftOrigin:    The rotation vector that contains the angels in degrees [ax+,ax-,ay+,ay-,az+,az-]
            interpolator:   The resampling filter interpolator. For gray images use sitk.sitkBSpline, and for binary images choose sitk.sitkNearestNeighbor
        Return:
            img_resampled:  The resampled image
    """
    
    T = sitk.Transform(3,sitk.sitkIdentity)

    resizeFilter = sitk.ResampleImageFilter()
    resizeFilter.SetTransform(T)

    oldSize      = img.GetSize()
    oldSpacing   = img.GetSpacing()
    
    newSize  = ( int(oldSize[0] * oldSpacing[0]  / newSpacing[0]),
        int(oldSize[1] * oldSpacing[1] / newSpacing[1]),
        int(oldSize[2] * oldSpacing [2] / newSpacing[2]) )

    oldOrigin    = img.GetOrigin()
    oldDirection = img.GetDirection()

    newOrigin    = [x + y for x, y in zip(oldOrigin, shiftOrigin)] 
    newDirection = oldDirection

    resizeFilter.SetOutputDirection(newDirection)
    resizeFilter.SetInterpolator(interpolator)
    resizeFilter.SetOutputSpacing(newSpacing)
    resizeFilter.SetOutputOrigin(newOrigin)
    resizeFilter.SetDefaultPixelValue(0)
    resizeFilter.SetSize(newSize)
    # resizeFilter.DebugOn()

    img_resampled = resizeFilter.Execute(img)
  
    return img_resampled

  
def transform_affine_sitk(fixed_image_sitk,translation_vector=[0,0,0],rotation_vector=[0,0,0],scaling_vector=[1,1,1],interpolator=sitk.sitkBSpline,spacing=None):
    """ This function transforms a nifti image
        Attributes:
            fixed_image_sitk:       The fixed image that will be transformed (simpleitk type)
            translation_vector:     The translation vector in mm [tx,ty,tz]
            rotation_angels:        The rotation vector that contains the angels in degrees [Rx,Ry,Rz]
            scaling_vector:         The scaling vector [Sx,Sy,Sz]
            interpolator:           The resampling filter interpolator. For gray images use sitk.sitkBSpline, and for binary images choose sitk.sitkNearestNeighbor
        Return:
            moving_image_sitk:      The moving image that has been transformed
    """

    # DefaultPixelValue = fixed_image_sitk.GetPixel(0,0,0)
    size        = fixed_image_sitk.GetSize()
    origin      = fixed_image_sitk.GetOrigin()
    direction   = fixed_image_sitk.GetDirection()

    if spacing is None:
        spacing = fixed_image_sitk.GetSpacing()
    
    dimension   = 3 
    Trans       = sitk.Transform(dimension, sitk.sitkAffine)

    # Translation
    tx, ty, tz = translation_vector/np.array(spacing)
    dt = np.array([tx,ty,tz,1])

    # Rotation
    theta_x, theta_y, theta_z = (mt.pi/180)*np.array(rotation_vector)
    
    Rx = np.array([
        [1, 0, 0],
        [0, mt.cos(theta_x), -mt.sin(theta_x)],
        [0, mt.sin(theta_x), mt.cos(theta_x)]])
    Ry = np.array([
        [mt.cos(theta_y), 0, mt.sin(theta_y)],
        [0, 1, 0],
        [-mt.sin(theta_y), 0, mt.cos(theta_y)]])
    Rz = np.array([
        [mt.cos(theta_z), -mt.sin(theta_z), 0],
        [mt.sin(theta_z), mt.cos(theta_z), 0],
        [0, 0, 1]])

    R = Rz.dot(Ry.dot(Rx.T).T)

    # Scale
    Sx, Sy, Sz = scaling_vector

    R[0,0] = R[0,0]/Sx
    R[1,1] = R[1,1]/Sy
    R[2,2] = R[2,2]/Sz
    
    # update transformation vector
    trans_vector = np.concatenate((R.flatten(),dt),axis=0)
    Trans.SetParameters(trans_vector)

    # print(Trans)

    # resample filter 
    resampleFilter = sitk.ResampleImageFilter()
    resampleFilter.SetTransform(Trans)
    resampleFilter.SetOutputDirection(direction)
    resampleFilter.SetInterpolator(interpolator)
    resampleFilter.SetOutputSpacing(spacing)
    resampleFilter.SetOutputOrigin(origin)
    resampleFilter.SetDefaultPixelValue(0)
    resampleFilter.SetSize(size)

    # transform the image
    moving_image_sitk = resampleFilter.Execute(fixed_image_sitk)

    return moving_image_sitk


def transform_skew_sitk(fixed_image_sitk, translation=(0,0,0), scale=(1,1,1), skew=(0,0,0,0,0,0), interpolator=sitk.sitkBSpline, spacing=None):
    """ This function transforms a nifti image
        Attributes:
            fixed_image_sitk:   The fixed image that will be transformed (simpleitk type)
            translation:        The translation vector in mm [tx,ty,tz]
            skew:               The rotation vector that contains the angels in degrees [ax+,ax-,ay+,ay-,az+,az-]
            scale:              The scaling vector [Sx,Sy,Sz]
            interpolator:       The resampling filter interpolator. For gray images use sitk.sitkBSpline, and for binary images choose sitk.sitkNearestNeighbor
        Return:
            moving_image_sitk:      The moving image that has been transformed
    """

    # DefaultPixelValue = fixed_image_sitk.GetPixel(0,0,0)
    size        = fixed_image_sitk.GetSize()
    origin      = fixed_image_sitk.GetOrigin()
    
    direction   = fixed_image_sitk.GetDirection()
    if spacing is None:
        spacing = fixed_image_sitk.GetSpacing()
    
    dimension   = 3 
    
    print skew
    print np.array(skew)
    skew        = np.tan(np.radians(np.array(skew))) #six eqaully spaced values in[0,1], an arbitrary choice
    versor      = (0,0,0,1.0)
    
    skewTransformer = sitk.ScaleSkewVersor3DTransform(scale, skew, versor, translation)

    print skewTransformer

    # resample filter 
    resampleFilter = sitk.ResampleImageFilter()
    resampleFilter.SetTransform(skewTransformer)
    resampleFilter.SetOutputDirection(direction)
    resampleFilter.SetInterpolator(interpolator)
    resampleFilter.SetOutputSpacing(spacing)
    resampleFilter.SetOutputOrigin(origin)
    resampleFilter.SetDefaultPixelValue(0)
    resampleFilter.SetSize(size)

    # transform the image
    moving_image_sitk = resampleFilter.Execute(fixed_image_sitk)

    return moving_image_sitk


def moveImages(fixed_image_sitk, shiftOrigin=(0,0,0), newSpacing=(1.,1.,1.), interpolator=sitk.sitkBSpline, translation_vector_1=[0,0,0], translation_vector_2=[0,0,0],
    skew_vector_1=(0,0,0,0,0,0), skew_vector_2=(0,0,0,0,0,0), scaling_vector_1=(1,1,1), scaling_vector_2=(1,1,1)):

    # reample image
    fixed_image_sitk0  = resample_sitk(fixed_image_sitk, newSpacing=newSpacing, shiftOrigin=shiftOrigin, interpolator=interpolator)
    fixed_image_sitk1  = resample_sitk(fixed_image_sitk, newSpacing=newSpacing, shiftOrigin=shiftOrigin, interpolator=interpolator)
    fixed_image_sitk2  = resample_sitk(fixed_image_sitk, newSpacing=newSpacing, shiftOrigin=shiftOrigin, interpolator=interpolator)
    
    # Transform image 
    moving_image_sitk0 = fixed_image_sitk0
    moving_image_sitk1 = transform_skew_sitk(fixed_image_sitk1, translation=translation_vector_1, scale=scaling_vector_1, skew=skew_vector_1, interpolator=interpolator)
    moving_image_sitk2 = transform_skew_sitk(fixed_image_sitk2, translation=translation_vector_2, scale=scaling_vector_2, skew=skew_vector_2, interpolator=interpolator)
    
    moving_image0      = sitk.GetArrayFromImage(moving_image_sitk0)    
    moving_image1      = sitk.GetArrayFromImage(moving_image_sitk1)
    moving_image2      = sitk.GetArrayFromImage(moving_image_sitk2)

    return moving_image_sitk0, moving_image0, moving_image1, moving_image2



def main():
    
    image_path  = 'brain_img.nii.gz'
    mask_path   = 'brain_mask.nii.gz'
    
    # read image
    fixed_image_sitk        = sitk.ReadImage(image_path)
    fixed_mask_sitk         = sitk.ReadImage(mask_path)
    fixed_image_mask_sitk   = sitk.Mask(fixed_image_sitk,fixed_mask_sitk)

    # reample image
    img_array = sitk.GetArrayFromImage(fixed_image_mask_sitk)

    threshold_filter = sitk.BinaryThresholdImageFilter()
    threshold_filter.SetLowerThreshold(1)
    threshold_filter.SetUpperThreshold((int)(img_array.max()))


    direction_ax = (
    1.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 
    0.0, 0.0, 1.0)

    direction_co = (
        1.0, 0.0, 0.0, 
        0.0, 0.0, 1.0, 
        0.0, 1.0, 0.0)

    direction_sa = (
        0.0, 0.0, 1.0, 
        0.0, 1.0, 0.0,
        1.0, 0.0, 0.0)


    img_ax      = img_array.copy()
    img_ax_sitk = sitk.GetImageFromArray(img_ax)
    img_ax_sitk.SetDirection(direction_ax)
    sitk.WriteImage(img_ax_sitk,'img_ax.nii.gz')
    img_ax_mask = threshold_filter.Execute(img_ax_sitk)
    sitk.WriteImage(img_ax_mask,'img_ax_mask.nii.gz')
    
    img_co      = np.swapaxes(img_array,0,1)
    img_co_sitk = sitk.GetImageFromArray(img_co)
    img_co_sitk.SetDirection(direction_co)
    sitk.WriteImage(img_co_sitk,'img_co.nii.gz')
    img_co_mask = threshold_filter.Execute(img_co_sitk)
    sitk.WriteImage(img_co_mask,'img_co_mask.nii.gz')

    img_sa      = np.swapaxes(img_array,0,2)
    img_sa_sitk = sitk.GetImageFromArray(img_sa)
    img_sa_sitk.SetDirection(direction_sa)
    sitk.WriteImage(img_sa_sitk,'img_sa.nii.gz')
    img_sa_mask = threshold_filter.Execute(img_sa_sitk)
    sitk.WriteImage(img_sa_mask,'img_sa_mask.nii.gz')


    # # randomize transformation parameters
    # translation_vector = np.random.choice(range(-4,4),3)      # random translations between [-4,4] mm
    # rotation_vector    = np.random.choice(range(-10,10),3)    # random rotation angels between [-10,10] degrees
    
    skew_unit_vector     = (1,1,1,1,1,1)

    # first transformation parameters
    skew_vector_1        = [1*x for x in skew_unit_vector]
    scaling_vector_1     = (1,1,1)
    translation_vector_1 = (0,0,0)   # translate x mm in x,y,z
    rotation_vector_1    = (0,0,0)   # rotate x degrees in x,y,z
    # first transformation parameters
    # skew_unit_vector     = (0,0,0,1,0,0)
    skew_vector_2        = [-1*x for x in skew_unit_vector]
    scaling_vector_2     = (1,1,1)
    translation_vector_2 = (0,0,0)    # translate x mm in x,y,z
    rotation_vector_2    = (0,0,0)    # rotate x degrees in x,y,z

    ################################################################################################################################################
    ## ---------------------------------------------------- Axial Image ------------------------------------------------------------------------- ##
    ################################################################################################################################################
    # axial-1 (a1) [0,1,2]
    newSpacing  = (1.5,1.5,3)
    shiftOrigin = (0,0,0)

    moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_ax_sitk, shiftOrigin=shiftOrigin, 
        newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
        skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_ax_mask, shiftOrigin=shiftOrigin, 
        newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
        skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    end = moving_image0.shape[0]

    sitk.WriteImage(moving_image_sitk0,'moving_img_a1.nii.gz')

    # sample slices -------------------------------------------------------------------------        
    moving_mask                 = moving_mask0
    moving_mask[range(0,end,3)] = moving_mask0[range(0,end,3)]
    moving_mask[range(1,end,3)] = moving_mask1[range(1,end,3)]
    moving_mask[range(2,end,3)] = moving_mask2[range(2,end,3)]
    moving_mask_sitk_final      = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_a1.nii.gz")

    moving_image                 = moving_image0    
    moving_image[range(0,end,3)] = moving_image0[range(0,end,3)]
    moving_image[range(1,end,3)] = moving_image1[range(1,end,3)]
    moving_image[range(2,end,3)] = moving_image2[range(2,end,3)]
    moving_image[moving_image<0] = 0

    moving_image_sitk_final      = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final      = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"img_a1.nii.gz")
    sitk.WriteImage(moving_image_sitk_final,"moving_image_a1.nii.gz")

    moving_image_sitk_final1     = sitk.GetImageFromArray(moving_image1)
    moving_image_sitk_final1.CopyInformation(moving_image_sitk0)
    sitk.WriteImage(moving_image_sitk_final1,"img_a2.nii.gz")

    moving_image_sitk_final2     = sitk.GetImageFromArray(moving_image2)
    moving_image_sitk_final2.CopyInformation(moving_image_sitk0)
    sitk.WriteImage(moving_image_sitk_final2,"img_a3.nii.gz")
    
    # # axial-2 (a2) [2,1,0]
    # shiftOrigin = (0,0,0)

    # moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_ax_sitk, shiftOrigin=shiftOrigin, 
    #     newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    # moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_ax_mask, shiftOrigin=shiftOrigin,
    #     newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    # end = moving_image0.shape[0]

    moving_mask                 = moving_mask0
    moving_mask[range(2,end,3)] = moving_mask0[range(2,end,3)]
    moving_mask[range(0,end,3)] = moving_mask1[range(0,end,3)]
    moving_mask[range(1,end,3)] = moving_mask2[range(1,end,3)]
    moving_mask_sitk_final      = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_a2.nii.gz")

    moving_image                 = moving_image0
    moving_image[range(2,end,3)] = moving_image0[range(2,end,3)]
    moving_image[range(0,end,3)] = moving_image1[range(0,end,3)]
    moving_image[range(1,end,3)] = moving_image2[range(1,end,3)]
    moving_image[moving_image<0] = 0
    moving_image_sitk_final      = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final      = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"moving_image_a2.nii.gz")

    # # axial-3 (a3) [1,2,0]
    # shiftOrigin = (0,0,0)

    # moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_ax_sitk, shiftOrigin=shiftOrigin,
    #     newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    # moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_ax_mask, shiftOrigin=shiftOrigin,
    #     newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    # end = moving_image0.shape[0]

    moving_mask                 = moving_mask0
    moving_mask[range(1,end,3)] = moving_mask0[range(1,end,3)]
    moving_mask[range(2,end,3)] = moving_mask1[range(2,end,3)]
    moving_mask[range(0,end,3)] = moving_mask2[range(0,end,3)]
    moving_mask_sitk_final      = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_a3.nii.gz")

    moving_image                 = moving_image0
    moving_image[range(1,end,3)] = moving_image0[range(1,end,3)]
    moving_image[range(2,end,3)] = moving_image1[range(2,end,3)]
    moving_image[range(0,end,3)] = moving_image2[range(0,end,3)]
    moving_image[moving_image<0] = 0
    moving_image_sitk_final      = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final     = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"moving_image_a3.nii.gz")


    ################################################################################################################################################
    ## ---------------------------------------------------- coronal Image ----------------------------------------------------------------------- ##
    ################################################################################################################################################
    newSpacing  = (1.5,3,1.5)
    shiftOrigin = (0,0,0)

    moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_co_sitk, shiftOrigin=shiftOrigin,
        newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
        skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_co_mask, shiftOrigin=shiftOrigin,
        newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
        skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    end = moving_image0.shape[0]
    
    sitk.WriteImage(moving_image_sitk0,'moving_img_c1.nii.gz')

    # sample slices -------------------------------------------------------------------------
    # coronal-1 (c1) [0,1,2]
    moving_mask                 = moving_mask0
    moving_mask[range(0,end,3)] = moving_mask0[range(0,end,3)]
    moving_mask[range(1,end,3)] = moving_mask1[range(1,end,3)]
    moving_mask[range(2,end,3)] = moving_mask2[range(2,end,3)]
    moving_mask_sitk_final      = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_c1.nii.gz")

    moving_image                 = moving_image0
    moving_image[range(0,end,3)] = moving_image0[range(0,end,3)]
    moving_image[range(1,end,3)] = moving_image1[range(1,end,3)]
    moving_image[range(2,end,3)] = moving_image2[range(2,end,3)]
    moving_image[moving_image<0] = 0
    moving_image_sitk_final      = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final      = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"moving_image_c1.nii.gz")

    # # coronal-2 (c2) [2,0,1]
    # shiftOrigin    = (0,0,0)

    # moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_co_sitk, shiftOrigin=shiftOrigin, 
    #     newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    # moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_co_mask, shiftOrigin=shiftOrigin, 
    #     newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    # end = moving_image0.shape[1]

    moving_mask                 = moving_mask0
    moving_mask[range(2,end,3)] = moving_mask0[range(2,end,3)]
    moving_mask[range(0,end,3)] = moving_mask1[range(0,end,3)]
    moving_mask[range(1,end,3)] = moving_mask2[range(1,end,3)]
    moving_mask_sitk_final      = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_c2.nii.gz")

    moving_image                 = moving_image0
    moving_image[range(2,end,3)] = moving_image0[range(2,end,3)]
    moving_image[range(0,end,3)] = moving_image1[range(0,end,3)]
    moving_image[range(1,end,3)] = moving_image2[range(1,end,3)]
    moving_image[moving_image<0] = 0
    moving_image_sitk_final      = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final      = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"moving_image_c2.nii.gz")

    # # coronal-3 (c3) [1,2,0]
    # shiftOrigin    = (0,0,0)

    # moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_co_sitk, shiftOrigin=shiftOrigin, 
    #     newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    # moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_co_mask, shiftOrigin=shiftOrigin, 
    #     newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    # end = moving_image0.shape[1]

    moving_mask                     = moving_mask0
    moving_mask[range(1,end,3)] = moving_mask0[range(1,end,3)]
    moving_mask[range(2,end,3)] = moving_mask1[range(2,end,3)]
    moving_mask[range(0,end,3)] = moving_mask2[range(0,end,3)]
    moving_mask_sitk_final          = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_c3.nii.gz")

    moving_image                    = moving_image0
    moving_image[range(1,end,3)] = moving_image0[range(1,end,3)]
    moving_image[range(2,end,3)] = moving_image1[range(2,end,3)]
    moving_image[range(0,end,3)] = moving_image2[range(0,end,3)]
    moving_image[moving_image<0]    = 0
    moving_image_sitk_final         = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final         = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"moving_image_c3.nii.gz")


    ################################################################################################################################################
    ## ---------------------------------------------------- saggital Image ---------------------------------------------------------------------- ##
    ################################################################################################################################################

    newSpacing  = (3,1.5,1.5)
    shiftOrigin = (0,0,0)

    moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_sa_sitk, shiftOrigin=shiftOrigin,
        newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
        skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_sa_mask, shiftOrigin=shiftOrigin,
        newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
        skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    end = moving_image0.shape[0]
    
    sitk.WriteImage(moving_image_sitk0,'moving_img_s1.nii.gz')

    # sample slices -------------------------------------------------------------------------
    # saggital-1 (s1) [0,1,2]
    moving_mask                 = moving_mask0
    moving_mask[range(0,end,3)] = moving_mask0[range(0,end,3)]
    moving_mask[range(1,end,3)] = moving_mask1[range(1,end,3)]
    moving_mask[range(2,end,3)] = moving_mask2[range(2,end,3)]
    moving_mask_sitk_final      = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_s1.nii.gz")

    moving_image                 = moving_image0
    moving_image[range(0,end,3)] = moving_image0[range(0,end,3)]
    moving_image[range(1,end,3)] = moving_image1[range(1,end,3)]
    moving_image[range(2,end,3)] = moving_image2[range(2,end,3)]
    moving_image[moving_image<0] = 0
    moving_image_sitk_final      = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final      = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"moving_image_s1.nii.gz")

    # # saggital-2 (s2) [2,0,1]
    # shiftOrigin = (0,0,0)

    # moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_sa_sitk, shiftOrigin=shiftOrigin,
    #     newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    # moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_sa_mask, shiftOrigin=shiftOrigin,
    #     newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    # end = moving_image0.shape[2]

    moving_mask                     = moving_mask0
    moving_mask[range(2,end,3)] = moving_mask0[range(2,end,3)]
    moving_mask[range(0,end,3)] = moving_mask1[range(0,end,3)]
    moving_mask[range(1,end,3)] = moving_mask2[range(1,end,3)]
    moving_mask_sitk_final          = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_s2.nii.gz")

    moving_image                    = moving_image0
    moving_image[range(2,end,3)] = moving_image0[range(2,end,3)]
    moving_image[range(0,end,3)] = moving_image1[range(0,end,3)]
    moving_image[range(1,end,3)] = moving_image2[range(1,end,3)]
    moving_image[moving_image<0]    = 0
    moving_image_sitk_final         = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final         = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"moving_image_s2.nii.gz")

    # # saggital-3 (s3) [1,2,0]
    # shiftOrigin = (0,0,0)

    # moving_image_sitk0, moving_image0, moving_image1, moving_image2 = moveImages(img_sa_sitk, shiftOrigin=shiftOrigin,
    #     newSpacing=newSpacing, interpolator=sitk.sitkBSpline, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)

    # moving_mask_sitk0, moving_mask0, moving_mask1, moving_mask2 = moveImages(img_sa_mask, shiftOrigin=shiftOrigin,
    #     newSpacing=newSpacing, interpolator=sitk.sitkNearestNeighbor, translation_vector_1=translation_vector_1, translation_vector_2=translation_vector_2,
    #     skew_vector_1=skew_vector_1, skew_vector_2=skew_vector_2, scaling_vector_1=scaling_vector_1, scaling_vector_2=scaling_vector_2)
    
    # end = moving_image0.shape[2]

    moving_mask                     = moving_mask0
    moving_mask[range(1,end,3)] = moving_mask0[range(1,end,3)]
    moving_mask[range(2,end,3)] = moving_mask1[range(2,end,3)]
    moving_mask[range(0,end,3)] = moving_mask2[range(0,end,3)]
    moving_mask_sitk_final          = sitk.GetImageFromArray(np.array(moving_mask,dtype=np.uint8))
    moving_mask_sitk_final.CopyInformation(moving_mask_sitk0)
    sitk.WriteImage(moving_mask_sitk_final,"moving_mask_s3.nii.gz")

    moving_image                    = moving_image0
    moving_image[range(1,end,3)] = moving_image0[range(1,end,3)]
    moving_image[range(2,end,3)] = moving_image1[range(3,end,3)]
    moving_image[range(0,end,3)] = moving_image2[range(0,end,3)]
    moving_image[moving_image<0]    = 0
    moving_image_sitk_final         = sitk.GetImageFromArray(moving_image)
    moving_image_sitk_final.CopyInformation(moving_image_sitk0)
    moving_image_sitk_final         = sitk.Mask(moving_image_sitk_final, moving_mask_sitk_final, 0)
    sitk.WriteImage(moving_image_sitk_final,"moving_image_s3.nii.gz")


if __name__ == '__main__':
    main()