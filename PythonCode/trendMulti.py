from PIL import Image, ImageDraw, ImageFont
import fitz  # PyMuPDF
import numpy as np

# -----------------------------
# Input files
# -----------------------------
top_path = "yearOutput.png"
bottom_pdf_path = "publishTrend.png"

# -----------------------------
# Helpers
# -----------------------------
def whiten_dark_background(img, threshold=40):
    arr = np.array(img.convert("RGB"))
    mask = (arr[:, :, 0] < threshold) & (arr[:, :, 1] < threshold) & (arr[:, :, 2] < threshold)
    arr[mask] = [255, 255, 255]
    return Image.fromarray(arr)

def pad_to_width(img, width, fill="white"):
    img = img.convert("RGB")
    if img.width == width:
        return img
    canvas = Image.new("RGB", (width, img.height), fill)
    x = (width - img.width) // 2
    canvas.paste(img, (x, 0))
    return canvas

def get_bold_font(size=90):
    font_candidates = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/Library/Fonts/Arial Bold.ttf",
        "C:/Windows/Fonts/arialbd.ttf",
    ]
    for fp in font_candidates:
        try:
            return ImageFont.truetype(fp, size)
        except Exception:
            continue
    return ImageFont.load_default()

# -----------------------------
# Load top PNG and clean background
# -----------------------------
top_img = Image.open(top_path).convert("RGB")
#top_img = whiten_dark_background(top_img)

# -----------------------------
# Render first page of PDF to image
# -----------------------------
bottom_img = Image.open(bottom_pdf_path).convert("RGB")
#bottom_img = whiten_dark_background(bottom_img)

# -----------------------------
# Pad both images to same width
# -----------------------------
target_width = max(top_img.width, bottom_img.width)
top_img = pad_to_width(top_img, target_width)
bottom_img = pad_to_width(bottom_img, target_width)

# -----------------------------
# Add white top margins for labels
# -----------------------------
label_margin = 120
top_with_margin = Image.new("RGB", (top_img.width, top_img.height + label_margin), "white")
top_with_margin.paste(top_img, (0, label_margin))

bottom_with_margin = Image.new("RGB", (bottom_img.width, bottom_img.height + label_margin), "white")
bottom_with_margin.paste(bottom_img, (0, label_margin))

# -----------------------------
# Combine vertically
# -----------------------------
gap = 40
combined = Image.new(
    "RGB",
    (target_width, top_with_margin.height + gap + bottom_with_margin.height),
    "white"
)
combined.paste(top_with_margin, (0, 0))
combined.paste(bottom_with_margin, (0, top_with_margin.height + gap))

# -----------------------------
# Remove black-bar artifact rows if any
# -----------------------------
#arr = np.array(combined)
#black_threshold = 40
#row_black_ratio = (
#    (arr[:, :, 0] < black_threshold) &
#    (arr[:, :, 1] < black_threshold) &
#    (arr[:, :, 2] < black_threshold)
#).mean(axis=1)

#artifact_rows = np.where(row_black_ratio > 0.5)[0]
#if len(artifact_rows) > 0:
#    mask = np.ones(arr.shape[0], dtype=bool)
#    mask[artifact_rows] = False
#    arr = arr[mask]
#    combined = Image.fromarray(arr)

# -----------------------------
# Add labels LAST so they are never removed
# -----------------------------
draw = ImageDraw.Draw(combined)
font = get_bold_font(size=90)

# Panel A label
draw.text((40, 15), "A", fill="black", font=font)

# Panel B label
b_y = top_with_margin.height + gap + 15
draw.text((40, b_y), "B", fill="black", font=font)

# -----------------------------
# Save final image
# -----------------------------
combined.save("trendMulti.png", dpi=(300, 300))
print("Saved as trendMulti.png")