import processing.core.PApplet;
import processing.core.PImage;
import processing.event.MouseEvent;

import java.util.Arrays;

public class Semiauto extends PApplet {
    String imgName = "/home/lev/IdeaProjects/allocator/modules/core/src/com/haulmont/arthur/allocator/core/solver_utils/tests/photoelast/data/Р01-125кг.jpg";

    int
            N = 4; // number of curve
    float
            sampleW = 52, // mm
            h = 5.7f,  // mm
            fSigma = 19.45659f / 10, // ??
            a = (34.745f - 16.388f) / 2f,
            b = sampleW / 2, // mm
            sigma = 125 /* / sampleW*/; // kg/mm^2

    PImage srcImg, modImg;

    boolean[] selectedPoints;
    int[] candidates;

    public void settings() {
        size(1000, 1000);
    }

    public void setup() {
        srcImg = loadImage(imgName);
        modImg = loadImage(imgName);

        srcImg.loadPixels();
        modImg.loadPixels();

        imgW = srcImg.width;
        imgH = srcImg.height;

        selectedPoints = new boolean[srcImg.width * srcImg.height];
        candidates = new int[srcImg.width * srcImg.height];

        //size(1000, 1000);

        estimations();
    }

    float imgX = 0, imgY = 0, imgW, imgH;

    int movingPoint = -1;

    public void draw() {
        background(0);
        image(modImg, imgX, imgY, imgW, imgH);

        fill(255);
        textSize(30);
        text(nfs(KI,0,3)
                + " min=" + nfs(min_KI,0,3)
                + " max=" + nfs(max_KI, 0, 3)
                + " abs=" + nfs(absDev_KI, 0, 3)
                + " std=" + nfs(std_KI, 0, 3), 10, 35);
        //text(KI, 10, 35);

//        fill(0, 0, 255);
//        text(KI_lastPoint, 10, 70);

/*
        if (centerX > 0) {
            fill(0, 0, 255);
            text(centerX, width - 300, 35);
        }
*/

        fill(255);
        for (int i = 0; i < F.length; i++) {
            text(sigma * sqrt(PI * a) * F[i], 10, 105 + i * 35);
        }

        text("SL="+softenLevel, 900, 35);
        text("SR="+softenRounds, 900, 70);
        text("MT="+methodsThreshold, 900, 105);
        text("EQ="+(equality ? "T" : "F"), 900, 140);
        text("N="+N, 900, 210);

        if (dragMode) {
            cursor();
        } else {
            noCursor();
            stroke(255);
            noFill();
            ellipse(mouseX, mouseY,
                    2 * paintR * imgW / srcImg.width,
                    2 * paintR * imgW / srcImg.width);
            fill(255);
        }
 /*
  if (movingPoint < 0) return;
  
  int 
    index = movingPoint - srcImg.width - 1, 
    maxIndex = srcImg.width * srcImg.height - 1;
  
  for (int i=0; i<9; i++) {
        
  }*/
    }

    boolean dragMode = true;
    int paintR = 10;

    private void selectPixel(int ofs) {
        if (candidates[ofs] <= methodsThreshold) return;
        selectedPoints[ofs] = true;
        modImg.pixels[ofs] = 0xFFFF0000;
    }

    private void paint(int xC, int yC) {
        int W = modImg.width;
        for (int y = - paintR; y <= 0; y++) {
            for (int x = - round(sqrt(paintR*paintR - y*y)); x <= 0; x++) {
                selectPixel((yC + y) * W + (xC + x));
                selectPixel((yC + y) * W + (xC - x));
                selectPixel((yC - y) * W + (xC + x));
                selectPixel((yC - y) * W + (xC - x));
            }
        }
        modImg.updatePixels();
    }

    public void mouseDragged() {
        if (dragMode) {
            imgX += mouseX - pmouseX;
            imgY += mouseY - pmouseY;
        } else if (mouseButton == LEFT) {
            int
                    px = round((mouseX - imgX) / imgW * srcImg.width),
                    py = round((mouseY - imgY) / imgH * srcImg.height);
            paint(px, py);
        }
    }

    public void mouseWheel(MouseEvent event) {
        if (dragMode) {
            float
                    oldW = imgW,
                    oldH = imgH,
                    zoom = event.getCount() < 0 ? 1.1f : 0.9f;

            imgW *= zoom;
            imgH *= zoom;

            imgX = mouseX + (imgX - mouseX) * imgW / oldW;
            imgY = mouseY + (imgY - mouseY) * imgH / oldH;
        } else {
            if (event.getCount() < 0) paintR++;
            else if (paintR > 1) paintR--;
        }
    }

    float centerX = 0, centerY = 0;
    int oldCenter = -1;

    public void mouseClicked() {
        int
                px = round((mouseX - imgX) / imgW * srcImg.width),
                py = round((mouseY - imgY) / imgH * srcImg.height);

        if (px < 0 || py < 0 || px >= srcImg.width || py >= srcImg.height) return;

        movingPoint = py * modImg.width + px;

        boolean center = mouseButton == RIGHT;

        if (center) {
            centerX = px * sampleW / srcImg.width;
            centerY = py * sampleW / srcImg.width;
            if (oldCenter >= 0) modImg.pixels[oldCenter] = srcImg.pixels[oldCenter];
            oldCenter = movingPoint;
            modImg.pixels[movingPoint] = 0xFF0000FF;
        } else {
            selectedPoints[movingPoint] = !selectedPoints[movingPoint];
            modImg.pixels[movingPoint] = selectedPoints[movingPoint] ? 0xFFFF0000 : srcImg.pixels[movingPoint];
        }

        // 0xFF0000FF - blue
        // 0xFFFF0000 - red
        // modImg.pixels[movingPoint] = center ? 0xFF0000FF : 0xFFFF0000;//center ? #0000FF: 0xFF0000;

        modImg.updatePixels();

        //if (oldCenter > 0 && points > 0) calculate();
    }

    public void clearSelection() {
        modImg = srcImg;
        srcImg = loadImage(imgName);
        srcImg.loadPixels();

        Arrays.fill(selectedPoints, false);
        Arrays.fill(candidates, 0);
    }

    public void keyPressed() {
        if (key == 'z' || key == 'Z' || key == 'я' || key == 'Я') {
            calculate();
        } else if (key == 'c' || key == 'C' || key == 'с' || key == 'С') {
            clearSelection();
        } else if (key == 's' || key == 'S' || key == 'ы' || key == 'Ы') {
            buildSkeletons();
        } else if (key == '-') {
            softenLevel--;
        } else if (key == '+') {
            softenLevel++;
        } else if (key == '[') {
            softenRounds--;
        } else if (key == ']') {
            softenRounds++;
        } else if (key == '\'') {
            methodsThreshold--;
        } else if (key == '\\') {
            methodsThreshold++;
        } else if (key == '/') {
            equality = !equality;
        } else if (key == ',' || key == 'б' || key == 'Б' || key == '<') {
            N--;
        } else if (key == '.' || key == '>' || key == 'ю' || key == 'Ю') {
            N++;
        } else if (key == 'd' || key == 'D') {
            dragMode = !dragMode;
        } else if (key == 'g' || key == 'G' || key == 'п' || key == 'П') {
            grayscale();
        }
    }

    float KI = 0, KI_lastPoint = 0, 
        min_KI = 0, max_KI = 0, absDev_KI = 0, std_KI = 0;

    void calculate() {
        KI = 0;

        min_KI = Float.MAX_VALUE;
        max_KI = -Float.MAX_VALUE;
        absDev_KI = 0;
        std_KI = 0;

        int points = 0;
        for (boolean selectedPoint : selectedPoints) {
            if (!selectedPoint) continue;
            points++;
        }

        float[] kiArr = new float[points];

        for (int i = 0, p = 0; i < selectedPoints.length; i++) {
            if (!selectedPoints[i]) continue;

            int
                    px = i % srcImg.width,
                    py = i / srcImg.width;

            float
                    pointX = px * sampleW / srcImg.width,
                    pointY = py * sampleW / srcImg.width,
                    dx = pointX - centerX,
                    dy = pointY - centerY,
                    r = sqrt(dx * dx + dy * dy),
                    sinPhi = -dy / r;

            float ki = sqrt(2 * PI * r) * N * fSigma / (h * sinPhi);
            KI += ki;
            kiArr[p++] = ki;

            if (ki > max_KI) max_KI = ki;
            if (ki < min_KI) min_KI = ki;
        }

        KI /= points;

        for (float v : kiArr) {
            absDev_KI += abs(v - KI);
            std_KI += sq(v - KI);
        }

        absDev_KI /= points;
        std_KI = sqrt(std_KI / points);
    }

    float[] F = new float[5];

    void estimations() {
        float ab = a / b;
        F[0] = sqrt(2 * b / (PI * a) * tan(PI * a / (2 * b)));
        F[1] = 1 + 0.128f * ab - 0.288f * ab * ab + 1.525f * ab * ab * ab;
        F[2] = sqrt(1.0f / cos(PI * ab / 2));
        F[3] = (1 - 0.5f * ab + 0.37f * ab * ab - 0.044f * ab * ab * ab) / sqrt(1 - ab);
        F[4] = (1 - 0.025f * ab * ab + 0.06f * ab * ab * ab * ab) * sqrt(1 / cos(PI * ab / 2));
        ;
    }

    public static void main(String[] passedArgs) {
        String[] appletArgs = new String[] { "Semiauto" };
        Semiauto semiauto = new Semiauto();
        PApplet.runSketch(appletArgs, semiauto);
    }

    public void soften(double[] p, int level) {
        if (level == 0) return;

        double[] r = p.clone();

        double window = level * 2d + 1d;

        for (int i = 0, maxI = p.length - 1; i <= maxI; i++) {
            double s = r[i];
            for (int j = 0; j < level; j++) {
                s += r[max(0,i - 1 - j)] + r[min(maxI, i + 1 + j)];
            }
            p[i] = s / window;
        }
    }

    int softenLevel = 4, softenRounds = 4, brightnessThreshold = 255*3 /*400*/;

    public static int pixelIntensity(int pixel) {
        return (pixel & 0xFF) + ((pixel >>> 8) & 0xFF) + ((pixel >>> 16) & 0xFF);
    }

    public void grayscale() {
        final int
                W = srcImg.width,
                H = srcImg.height;

        for (int i = 0; i < W*H; i++) {
            int pixel = srcImg.pixels[i];
            int p = (int)Math.round((pixelIntensity(pixel)) / 3d);
            srcImg.pixels[i] = p + (p << 8) + (p << 16) + 0xFF000000;
            modImg.pixels[i] = p + (p << 8) + (p << 16) + 0xFF000000;
        }

        srcImg.updatePixels();
        modImg.updatePixels();
    }

    public void buildSkeletons() {
        buildSkeletons_HorVer(candidates);
        buildSkeletons_Diagonal(candidates);
        buildSkeletons_RameshHorVer(candidates);
        buildSkeletons_RameshDiagonal(candidates);

        showResult(candidates);
    }

    int methodsThreshold = 3;

    public void showResult(int[] result) {
        for (int i = 0; i < result.length; i++) {
            if (result[i] > methodsThreshold) modImg.pixels[i] = 0xFFFFFFFF;
        }
        modImg.updatePixels();
    }

    boolean equality = false;

    public void buildSkeletons_HorVer(int[] result) {
        final int
                W = srcImg.width,
                H = srcImg.height;

        double[] p = new double[W];

        for (int y = 0, ptr = 0; y < H; y++) {
            int ptrRes = ptr;
            for (int x = 0; x < W; x++) {
                int pixel = srcImg.pixels[ptr++];
                p[x] = pixelIntensity(pixel);
            }

            for (int i = 0; i < softenRounds; i++) soften(p, softenLevel);

            for (int x = 1; x < W - 1; x++) {
                ptrRes++;
                if (equality) {
                    if (!(p[x] <= p[x - 1] && p[x] <= p[x + 1] /* && p[x] < brightnessThreshold */)) continue;
                } else {
                    if (!(p[x] < p[x - 1] && p[x] < p[x + 1] /* && p[x] < brightnessThreshold */)) continue;
                }
                result[ptrRes]++;
            }
        }

        p = new double[H];

        for (int x = 0; x < W; x++) {
            int ptr = x;
            for (int y = 0; y < H; y++) {
                int pixel = srcImg.pixels[ptr];
                ptr += W;
                p[y] = pixelIntensity(pixel);
            }

            for (int i = 0; i < softenRounds; i++) soften(p, softenLevel);

            for (int y = 1; y < H - 1; y++) {
                if (equality) {
                    if (!(p[y] <= p[y - 1] && p[y] <= p[y + 1] /* && p[x] < brightnessThreshold */)) continue;
                } else {
                    if (!(p[y] < p[y - 1] && p[y] < p[y + 1] /* && p[x] < brightnessThreshold */)) continue;
                }
                result[y * W + x]++;
            }
        }
    }

    public void buildSkeletons_Diagonal(int[] result) {
        final int
                W = srcImg.width,
                H = srcImg.height;

        for (int i = 0, x = 0, y = 0; i < W+H-1; i++) {

            int ptr = y * W + x;

            double[] p = new double[Math.min(y+1, W - x)];

            for (int j = 0; j < p.length; j++) {
                int pixel = srcImg.pixels[ptr];
                ptr -= W-1;
                p[j] = pixelIntensity(pixel);
            }

            for (int j = 0; j < softenRounds; j++) soften(p, softenLevel);

            for (int j = 1, xx = x, yy = y; j < p.length - 1; j++) {
                if (equality) {
                    if (p[j] <= p[j - 1] && p[j] <= p[j + 1]) result[yy * W + xx]++;
                } else {
                    if (p[j] < p[j - 1] && p[j] < p[j + 1]) result[yy * W + xx]++;
                }
                xx++;
                yy--;
            }

            if (y < H - 1) y++; else x++;
        }

        for (int i = 0, x = W - 1, y = 0; i < W+H-1; i++) {

            int ptr = y * W + x;

            double[] p = new double[Math.min(y+1, x+1)];

            for (int j = 0; j < p.length; j++) {
                int pixel = srcImg.pixels[ptr];
                ptr -= W+1;
                p[j] = pixelIntensity(pixel);
            }

            for (int j = 0; j < softenRounds; j++) soften(p, softenLevel);

            for (int j = 1, xx = x, yy = y; j < p.length - 1; j++) {
                if (equality) {
                    if (p[j] <= p[j - 1] && p[j] <= p[j + 1]) result[yy * W + xx]++;
                } else {
                    if (p[j] < p[j - 1] && p[j] < p[j + 1]) result[yy * W + xx]++;
                }
                xx--;
                yy--;
            }

            if (y < H - 1) y++; else x--;
        }
    }

    public void buildSkeletons_RameshHorVer(int[] result) {
        final int
                W = srcImg.width,
                H = srcImg.height;

        for (int y = 2; y < H-2; y++) {
            for (int x = 2; x < W-2; x++) {
                int ptr = y * W + x - 2;

                int leftCol =
                        pixelIntensity(srcImg.pixels[ptr - W]) +
                        pixelIntensity(srcImg.pixels[ptr]) +
                        pixelIntensity(srcImg.pixels[ptr + W]);

                ptr += 4;

                int rightCol =
                        pixelIntensity(srcImg.pixels[ptr - W]) +
                        pixelIntensity(srcImg.pixels[ptr]) +
                        pixelIntensity(srcImg.pixels[ptr + W]);

                ptr -= 2;

                int centerCol =
                        pixelIntensity(srcImg.pixels[ptr - W]) +
                        pixelIntensity(srcImg.pixels[ptr]) +
                        pixelIntensity(srcImg.pixels[ptr + W]);

                ptr -= W << 1;

                int topRow =
                        pixelIntensity(srcImg.pixels[ptr - 1]) +
                        pixelIntensity(srcImg.pixels[ptr]) +
                        pixelIntensity(srcImg.pixels[ptr + 1]);

                ptr += W << 2;

                int bottomRow =
                        pixelIntensity(srcImg.pixels[ptr - 1]) +
                        pixelIntensity(srcImg.pixels[ptr]) +
                        pixelIntensity(srcImg.pixels[ptr + 1]);

                ptr -= W << 1;

                int centerRow =
                        pixelIntensity(srcImg.pixels[ptr - 1]) +
                        pixelIntensity(srcImg.pixels[ptr]) +
                        pixelIntensity(srcImg.pixels[ptr + 1]);

                if (centerCol < leftCol && centerCol < rightCol) result[ptr]++;
                if (centerRow < topRow && centerRow < bottomRow) result[ptr]++;
            }
        }
    }

    public void buildSkeletons_RameshDiagonal(int[] result) {
        final int
                W = srcImg.width,
                H = srcImg.height;

        for (int y = 2; y < H-2; y++) {
            for (int x = 2; x < W-2; x++) {
                int ptr = y * W + x;

                double centerTL = (
                        pixelIntensity(srcImg.pixels[ptr - W - 1]) +
                        pixelIntensity(srcImg.pixels[ptr]) +
                        pixelIntensity(srcImg.pixels[ptr + W + 1])
                ) / 3d;

                double BL = (
                        pixelIntensity(srcImg.pixels[ptr + W - 2]) +
                        pixelIntensity(srcImg.pixels[ptr + 2*W - 1])
                ) / 2d;

                double TR = (
                        pixelIntensity(srcImg.pixels[ptr - W + 2]) +
                        pixelIntensity(srcImg.pixels[ptr - 2*W + 1])
                ) / 2d;

                double centerTR = (
                        pixelIntensity(srcImg.pixels[ptr - W + 1]) +
                        pixelIntensity(srcImg.pixels[ptr]) +
                        pixelIntensity(srcImg.pixels[ptr + W - 1])
                ) / 3d;

                double BR = (
                        pixelIntensity(srcImg.pixels[ptr + W + 2]) +
                        pixelIntensity(srcImg.pixels[ptr + 2*W + 1])
                ) / 2d;

                double TL = (
                        pixelIntensity(srcImg.pixels[ptr - W - 2]) +
                        pixelIntensity(srcImg.pixels[ptr - 2*W - 1])
                ) / 2d;

                if (centerTL < BL && centerTL < TR) result[ptr]++;
                if (centerTR < BR && centerTR < TL) result[ptr]++;
            }
        }
    }

}