/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package browser.table;

import browser.ICONS_MANAGER;
import browser.imageitems.TableImageItem;
import ij.ImagePlus;
import java.awt.Color;
import java.awt.Component;
import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.table.TableCellRenderer;

/**
 *
 * @author Juanjo Vega
 */
public class ImagesTableRenderer extends JLabel implements TableCellRenderer {

    protected final static Border BORDER_SELECTED = BorderFactory.createLineBorder(Color.RED, 1);
    protected final static Border BORDER_FOCUSED = BorderFactory.createLineBorder(Color.RED, 3);
    protected final static Color BACKGROUND = UIManager.getColor("Table.background");
    protected final static Color BACKGROUND_SELECTED = UIManager.getColor("Table.selectionBackground");
    protected final static Color FOREGROUND = UIManager.getColor("Table.foreground");
    protected final static Color FOREGROUND_SELECTED = UIManager.getColor("Table.selectionForeground");

    public ImagesTableRenderer() {
        super();
    }

    public Component getTableCellRendererComponent(JTable table, Object object, boolean isSelected, boolean hasFocus, int row, int column) {
        TableImageItem item = (TableImageItem) object;

        if (item != null) {
//            System.out.println(" *** Rendering: " + item.getLabel());

            // Loads image...
            ImagePlus img = item.getPreview();

            // ... and sets it.
            if (img != null) {
                setIcon(new ImageIcon(img.getImage()));
                img.close();
            } else {
                setIcon(ICONS_MANAGER.MISSING_ITEM);
            }

            setOpaque(true);
            setHorizontalAlignment(CENTER);
            setHorizontalTextPosition(CENTER);
            setVerticalTextPosition(BOTTOM);

            // Tooltip.
            setToolTipText(item.getLabel());

            // (Shows label only when required).
            if (((ImagesTableModel) table.getModel()).isShowingLabels()) {
                setText(item.getLabel());
            } else {
                setText(null);
            }

            if (item.isSelected()) {    // If is selected...
                if (hasFocus) {    // ...and focused as well.
                    setBorder(BORDER_FOCUSED);
                    setBackground(BACKGROUND_SELECTED);
                    setForeground(FOREGROUND_SELECTED);
                } else {    // ...otherwise.
                    setBorder(BORDER_SELECTED);
                    setBackground(BACKGROUND_SELECTED);
                    setForeground(FOREGROUND_SELECTED);
                }
            } else {
                setBorder(null);
                setBackground(BACKGROUND);
                setForeground(FOREGROUND);
            }
        } else {
            setIcon(null);
            setText(null);
            setBorder(null);
            setBackground(BACKGROUND);
            setForeground(FOREGROUND_SELECTED);
        }

        return this;
    }
}
