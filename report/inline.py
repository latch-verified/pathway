from pathlib import Path
import re, os, glob

if __name__ == "__main__":
    print("Inlining js and css files into html file...")

    htmlpath = Path("./build/index.html")

    with htmlpath.open("r") as fhtml:
        html = fhtml.read()

        with open(glob.glob("./build/static/css/*.css")[0], "r") as fcss:
            stylenode = f"<style>{fcss.read()}</style>"
            html = re.sub("<link.*?main.*?css.*?>", stylenode, html)

        with open(glob.glob("./build/static/js/*.js")[0], "r") as fjs:
            # Using type="module" to defer the script execution until after the
            # page loads. The defer keyword only works for scripts with a valid,
            # existing src attribute
            scriptnode = f'<script type="module">{fjs.read()}</script>'
            # Using lambda to stop Regex module from attempting to escape the
            # escape codes typically found in the JavaScript file
            html = re.sub("<script.*?main.*?js.*?</script>", lambda _: scriptnode, html)

    # Remove non-index.html build outputs
    os.system("rm -rf ./build/static/")
    os.remove("./build/asset-manifest.json")

    with htmlpath.open("w") as fhtml:
        fhtml.write(html)

    print("Completed inlining")
