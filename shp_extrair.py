import fiona
from shapely.geometry import shape, box

# Defina a área de interesse (bounding box) - você pode ajustar essas coordenadas

p1 = (-43.171140, -22.954110)
p2 = (-43.169465, -22.953343)

def extrarir_do_rio(p1,p2):
    minx, miny = min(p1[0], p2[0])-1/360, min(p1[1], p2[1])-1/360
    maxx, maxy = max(p1[0], p2[0])+1/360, max(p1[1], p2[1])+1/360
    bbox = box(minx, miny, maxx, maxy)

    nome = "construcoes"

    # Abra o shapefile grande
    input_shapefile = 'shapefiles\\'+nome+'.shp'
    output_shapefile = 'shapefiles\derivados\\'+nome+'.shp'

    # Crie um novo shapefile para armazenar os dados extraídos
    with fiona.open(input_shapefile, 'r') as src:
        meta = src.meta
        with fiona.open(output_shapefile, 'w', **meta) as dst:
            for feature in src:
                if feature['geometry'] is not None:
                    geom = shape(feature['geometry'])
                    if geom.intersects(bbox):
                        # Se a geometria intersecta a bounding box, adicioná-la ao novo shapefile
                        dst.write(feature)
    return output_shapefile