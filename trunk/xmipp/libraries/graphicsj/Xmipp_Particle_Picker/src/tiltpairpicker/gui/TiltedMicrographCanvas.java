package tiltpairpicker.gui;

import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.SwingUtilities;

import picker.gui.WindowUtils;
import picker.model.Particle;
import tiltpairpicker.model.ParticlePairPicker;
import tiltpairpicker.model.TiltedMicrograph;


public class TiltedMicrographCanvas extends ImageCanvas implements MouseListener{

	private ParticlePairPickerJFrame frame;
	private TiltedMicrograph tiltedmicrograph;
	private Particle dragged;
	private ParticlePairPicker pppicker;
	private ImageWindow iw;
	

	public TiltedMicrographCanvas(ParticlePairPickerJFrame frame) {
		super(frame.getUntiltedMicrograph().getTiltedMicrograph().getImage());
		this.tiltedmicrograph = frame.getUntiltedMicrograph().getTiltedMicrograph();
		this.frame = frame;
		this.pppicker = frame.getParticlePairPicker();
		iw = new ImageWindow(imp, this);
		WindowUtils.centerScreen(0.7, iw);
		
	}
	
	public void updateMicrograph() {
		this.tiltedmicrograph = frame.getUntiltedMicrograph().getTiltedMicrograph();
		imp = tiltedmicrograph.getImage();
		iw.setImage(imp);
		iw.setTitle(tiltedmicrograph.getName());
		setImageUpdated();
		repaint();
	}
	
	
	public void mouseEntered(MouseEvent e) {
		if (frame.getTool() != Tool.PICKER) {
			super.mouseEntered(e);
			return;
		}
		setCursor(crosshairCursor);
	}
	
	public void mouseMoved(MouseEvent e) {
		if (frame.getTool() != Tool.PICKER) {
			super.mouseMoved(e);
			return;
		}
		setCursor(crosshairCursor);
	}

	/**
	 * Adds particle or updates its position if onpick. If ondeletepick removes
	 * particle. Considers owner for selection to the first particle containing
	 * point. Sets dragged if onpick
	 */

	public void mousePressed(int x, int y) {
		setupScroll(x, y);
		
	}

	
	/**
	 * Updates particle position and repaints. Sets dragged to null at the end
	 */
	public void mouseReleased(MouseEvent e) {
		if (frame.getTool() != Tool.PICKER) {
			super.mouseReleased(e);
			return;
		}
		dragged = null;
	}

	/**
	 * Updates particle position and repaints if onpick.
	 */
	public void mouseDragged(int x, int y) {
			scroll(x, y);
		
	}
	
	public void mouseWheelMoved(int x, int y, int rotation) {
		
		if (rotation < 0)
			zoomIn(x, y);
		else
			zoomOut(x, y);

	}

	public void paint(Graphics g) {
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		
	}
	

	private void drawShape(Graphics2D g2, Particle p, int x0, int y0, int radius, boolean all) {
		
		int x = (int) ((p.getX() - x0) * magnification);
		int y = (int) ((p.getY() - y0) * magnification);
		int distance = (int)(10 * magnification);
		if (frame.isShapeSelected(Shape.Rectangle) || all)
			g2.drawRect(x - radius, y - radius, radius * 2, radius * 2);
		if (frame.isShapeSelected(Shape.Circle) || all)
			g2.drawOval(x - radius, y - radius, radius * 2, radius * 2);
		if (frame.isShapeSelected(Shape.Center) || all) 
		{
			g2.drawLine(x, y - distance, x, y + distance);
			g2.drawLine(x + distance, y, x - distance, y);
		}
	}

	

}
