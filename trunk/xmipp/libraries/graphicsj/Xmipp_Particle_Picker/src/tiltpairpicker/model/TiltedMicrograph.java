package tiltpairpicker.model;

import ij.ImagePlus;

import java.io.File;
import java.util.List;

import javax.swing.Icon;

import picker.model.Constants;
import picker.model.Family;
import picker.model.FamilyState;
import picker.model.MicrographFamilyData;
import picker.model.MicrographFamilyState;
import picker.model.Particle;

public class TiltedMicrograph {
	
	private String imagefile;
	private String name;
	private ImagePlus image;
	
	public TiltedMicrograph(String image) {
		
		this.imagefile = image;
		if(!new File(image).exists())
			throw new IllegalArgumentException(Constants.getNoSuchFieldValueMsg("image", image));
		this.name = getName(image);

	}

	
	public String getImageFile() {
		return imagefile;
	}


	public String getName() {
		return name;
	}
	
	public static String getName(String file)
	{
		String[] tokens = file.split(File.separator);
		return  tokens[tokens.length - 1];
	}
	
	

	public boolean isEmpty() {
		// TODO Auto-generated method stub
		return false;
	}

	
	public ImagePlus getImage()
	{
		if(image == null)
			image = new ImagePlus(imagefile);
		return image;
	}
	
	public void releaseImage()
	{
		image = null;
	}

}
