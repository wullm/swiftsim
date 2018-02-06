htmlfile = open("index.html", "w")

htmlfile.write("<!DOCTYPE html>")
htmlfile.write("<html lang=\"en\">")
htmlfile.write("<head><title>1D Sod shock convergence</title></head>")
htmlfile.write("<body><h1>1D Sod shock convergence</h1>")

schemes = {"gizmo": "GIZMO",
           "gadget2": "Gadget2",
           "hopkins": "Pressure-entropy SPH"
          }

ncell = [100, 200, 400, 800, 1600, 3200]

htmlfile.write("<p><img src=\"SodShock_1D_convergence.png\""
               " height=\"400px\" />")
htmlfile.write("<img src=\"SodShock_1D_timings.png\" height=\"400px\" /></p>")

htmlfile.write("<table>")
htmlfile.write("<tr>")
for scheme in schemes:
  htmlfile.write("<th>{0}</th>".format(schemes[scheme]))
htmlfile.write("</tr>")

for n in ncell:
  htmlfile.write("<tr>")
  for scheme in schemes:
    htmlfile.write("<td width=\"300px\">")
    htmlfile.write("<img src=\"result_{0}_{1}.png\" width=\"300px\" />".format(
                     scheme, n))
    htmlfile.write("</td>")
  htmlfile.write("</tr>")

htmlfile.write("</table>")
htmlfile.write("</body>")
htmlfile.write("</html>")
