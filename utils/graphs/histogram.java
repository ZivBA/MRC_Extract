package utils.graphs;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.Map;
import java.util.TreeMap;

/**
 * Created by zivben on 09/08/15.
 */
public class histogram {

	public histogram(double[][] data) {
		// For this example, I just randomised some data, you would
		// Need to load it yourself...
		//		int width = 256;
		//		int height = 256;

		//		int[][] data = new int[width][height];
		//		for (int c = 0; c < height; c++) {
		//			for (int r = 0; r < width; r++) {
		//				data[c][r] = (int) (256 * Math.random());
		//			}
		//		}


		Map<Double, Double> mapHistory = new TreeMap<Double, Double>();
		for (int c = 0; c < data.length; c++) {
			for (int r = 0; r < data[c].length; r++) {
				double value = data[c][r];
				Double amount;
				if (mapHistory.containsKey(value)) {
					amount = mapHistory.get(value);
					amount++;
				} else {
					amount = 0d;
				}
				mapHistory.put(value, amount);
			}
		}
		JFrame frame = new JFrame("Test");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setLayout(new BorderLayout());
		frame.add(new JScrollPane(new Graph(mapHistory)));
		frame.pack();
		frame.setLocationRelativeTo(null);
		frame.setVisible(true);
	}

	protected class Graph extends JPanel {

		protected static final int MIN_BAR_WIDTH = 4;
		private Map<Double, Double> mapHistory;

		public Graph(Map<Double, Double> mapHistory) {
			this.mapHistory = mapHistory;
			int width = (mapHistory.size() * MIN_BAR_WIDTH) + 11;
			Dimension minSize = new Dimension(width, 128);
			Dimension prefSize = new Dimension(width, 256);
			setMinimumSize(minSize);
			setPreferredSize(prefSize);
		}

		@Override
		protected void paintComponent(Graphics g) {
			super.paintComponent(g);
			if (mapHistory != null) {
				int xOffset = 5;
				int yOffset = 5;
				int width = getWidth() - 1 - (xOffset * 2);
				int height = getHeight() - 1 - (yOffset * 2);
				Graphics2D g2d = (Graphics2D) g.create();
				g2d.setColor(Color.DARK_GRAY);
				g2d.drawRect(xOffset, yOffset, width, height);
				int barWidth = Math.max(MIN_BAR_WIDTH,
						(int) Math.floor((float) width
								/ (float) mapHistory.size()));
				System.out.println("width = " + width + "; size = "
						+ mapHistory.size() + "; barWidth = " + barWidth);
				double maxValue = 0;
				for (Double key : mapHistory.keySet()) {
					Double value = mapHistory.get(key);
					maxValue = Math.max(maxValue, value);
				}
				int xPos = xOffset;
				for (Double key : mapHistory.keySet()) {
					Double value = mapHistory.get(key);
					long barHeight = Math.round((value
							/ maxValue) * height);
					g2d.setColor(new Color(key.intValue(), key.intValue(), key.intValue()));
					long yPos = height + yOffset - barHeight;
					//Rectangle bar = new Rectangle(xPos, yPos, barWidth, barHeight);
					Rectangle2D bar = new Rectangle2D.Float(
							xPos, yPos, barWidth, barHeight);
					g2d.fill(bar);
					g2d.setColor(Color.DARK_GRAY);
					g2d.draw(bar);
					xPos += barWidth;
				}
				g2d.dispose();
			}
		}
	}
}
