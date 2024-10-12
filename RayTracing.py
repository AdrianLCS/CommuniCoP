import geopandas as gpd
import mpmath
from shapely.geometry import LineString, Point, Polygon, MultiLineString, MultiPolygon, GeometryCollection
import numpy as np
import cmath
import Modelos
#from scipy.integrate import quad
from mpmath import quad, sin, pi, cos

a = 6378137  # m
b = 6356752  # m
mi0 = 4 * np.pi * 1e-7
e0 = 8.85e-12

azes = [+1.595769140,-0.000001702,-6.808568854,-0.000576361,+6.920691902,-0.016898657,-3.050485660,-0.075752419,+0.850663781,-0.025639041,-0.150230960,+0.034404779]
bzes = [-0.000000033,+4.255387524,-0.000092810,-7.780020400,-0.009520895,+5.075161298,-0.138341947,-1.363729124,-0.403349276,+0.702222016,-0.216195929,+0.019547031]
czes = [+0.000000000,-0.024933975,+0.000003936,+0.005770956,+0.000689892,-0.009497136,+0.011948809,-0.006748873,+0.000246420,+0.002102967,-0.001217930,+0.000233939]
dzes = [+0.199471140,+0.000000023,-0.009351341,+0.000023006,+0.004851466,+0.001903218,-0.017122914,+0.029064067,-0.027928955,+0.016497308,-0.005598515,+0.000838386]

def Fv(v):
    x = 0.5*np.pi*(v**2)
    if x < 4:
        var = 0
        for i in range(len(azes)):
            var = var+((azes[i]-1j*bzes[i])*((x/4)**i))
        return np.exp(1j*x)*((x/4)**0.5)*var
    else:
        var = 0
        for i in range(len(azes)):
            var=var+((czes[i]-1j*dzes[i])*((4/x)**i))
        return ((1+1j)/2)+np.exp(1j*x)*((4/x)**0.5)*var


def R(lat):
    """Essa função retorna o raio da terra em função da latitude. Usada para o cáculuo de distância entre dois pontos"""
    return (((((a ** 2) * np.cos(lat * np.pi / 180)) ** 2) + (((b ** 2) * np.sin(lat * np.pi / 180)) ** 2)) / (
            ((a * np.cos(lat * np.pi / 180)) ** 2) + ((b * np.sin(lat * np.pi / 180)) ** 2))) ** 0.5


def deg2rad(degrees):
    """Converte graus para radianos"""
    radians = degrees * np.pi / 180
    return radians


def getDistanceBetweenPointsNew(latitude1, longitude1, latitude2, longitude2):
    """Formula para cálculo de distância entre duas coordenadas"""
    theta = longitude1 - longitude2
    latitude1 = deg2rad(latitude1)
    latitude2 = deg2rad(latitude2)
    longitude1 = deg2rad(longitude1)
    longitude2 = deg2rad(longitude2)
    distance = 2 * R((latitude1 + latitude2) / 2) * np.arcsin(((np.sin((latitude2 - latitude1) / 2)) ** 2 +
                                                               np.cos(latitude1) * np.cos(latitude2) * ((np.sin(
                (longitude2 - longitude1) / 2)) ** 2)) ** 0.5)

    return distance


def load_shapefile(filepath):
    """Load shapefile and return GeoDataFrame."""
    return gpd.read_file(filepath)


def generate_ray_directions_difrac(direcao0, num_azimuths):  # direcao0 é o vetor uintário (p2-p1)/norm(p2-p1)
    """Generate equally spaced ray directions in 2D."""
    directions = []
    ang_dir0 = np.arccos(direcao0[0])
    azimuth_angles = np.linspace(ang_dir0 - np.pi / 4, ang_dir0 + np.pi / 4, num_azimuths, endpoint=False)

    for azimuth in azimuth_angles:
        x = np.cos(azimuth)
        y = np.sin(azimuth)
        directions.append(np.array([x, y]))  # Z-component is zero for 2D
    return directions


def generate_ray_directions(num_azimuths):
    """Generate equally spaced ray directions in 2D."""
    directions = []
    azimuth_angles = np.linspace(0, 2 * np.pi, num_azimuths, endpoint=False)

    for azimuth in azimuth_angles:
        x = np.cos(azimuth)
        y = np.sin(azimuth)
        directions.append((x, y))  # Z-component is zero for 2D
    return directions


def reflec_raio(segment, ray_position, itersec_point):
    print("""Calculate the normal of a line segment and orient it towards the ray position.""")
    itersec_point = np.array(itersec_point.x, itersec_point.y)

    if isinstance(segment, LineString):
        pontos = list(segment.coords)
        p1, p2 = np.array(pontos[0]), np.array(pontos[1])
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]

        raio = np.array(itersec_point - ray_position) / np.linalg.norm(np.array(itersec_point - ray_position))

        dir_parede = np.array([dx, dy])

        dir_parede = dir_parede / np.linalg.norm(dir_parede)
        normal = np.array([-dy, dx])
        normal = normal / np.linalg.norm(normal)

        mag_comp1 = np.dot(dir_parede, raio)
        mag_comp2 = np.dot(normal, raio)

        if mag_comp1 < 0:
            mag_comp1 = -mag_comp1
            dir_parede = -dir_parede

        # Adjust normal to point towards the ray source
        if mag_comp2 > 0:
            normal = -normal
        else:
            mag_comp2 = -mag_comp2

        raio_refletido = mag_comp1 * dir_parede + mag_comp2 * normal

        print(dir_parede)
        print(raio)
        print(raio_refletido)

        return normal, mag_comp2  # mag comp é o cosseno
    return None


def trace_rays(shapeData, tx_position, rx_position, num_azimuths=1000, max_reflections=3, max_diffractions=1,
               diffraction_distance=0.5):
    """Trace rays from tx_position to rx_position in 2D."""
    quinas = []
    ray_paths = []
    directions = generate_ray_directions(num_azimuths)
    shape_ja_testados = []
    tx_position = np.array(tx_position[:2])  # Only use XY components
    distancia = np.linalg.norm(rx_position[:2] - tx_position)
    max_distance = distancia * 2 * np.pi / num_azimuths #distancia * 2 * np.pi / num_azimuths
    dir_tx_rx = (rx_position[:2] - tx_position) / distancia
    directions.append(dir_tx_rx)

    cont0 = 0
    for direction in directions:
        ray_position = np.array(tx_position)
        path = []
        reflections = 0
        diffractions = 0
        hit_rx = False
        steps = 0  # Counter for number of steps
        cont0 += 1
        intersection_line = 0
        incidence_angle = 0  # cosseno do angulo do raio refletido com a normal
        print('raio: ' + str(cont0))

        cont1 = 0
        print('shapeData')

        while (reflections <= max_reflections) and not hit_rx and steps < 1:
            cont1 += 1
            print('subraio: ' + str(cont1))
            ray_direction = np.array(direction)
            ray_line = LineString([ray_position, ray_position + ray_direction * 0.1])
            intersections = []
            intersected_segments = []
            intersection_point = Point([0, 0])

            for shape in shapeData.geometry:
                intersecao = ray_line.intersection(shape)
                if not intersecao.is_empty:
                    intersected_segments.append(shape)

            if intersected_segments:
                closest_shape = min(intersected_segments, key=lambda seg: seg.distance(Point(ray_position)))

                intersection = ray_line.intersection(closest_shape)
                if intersection.geom_type == 'Point':
                    # print('interseção ponto ocorreu')
                    intersection_point = intersection
                    intersections.append(
                        (intersection, closest_shape))  # (ponto, reta) shape é uma linestrig nesse caso
                elif intersection.geom_type == 'LineString':
                    P1 = Point(intersection.coords[0][:2])
                    P2 = Point(intersection.coords[-1][:2])
                    if P1.distance(Point(ray_position)) < P2.distance(Point(ray_position)):
                        intersection_point = P1
                    else:
                        intersection_point = P2
                elif intersection.geom_type == 'MultiLineString':
                    Pn = []
                    Pn_dist = []
                    for string in intersection.geoms:
                        P1 = Point(string.coords[0][:2])
                        Pn.append(P1)
                        Pn_dist.append(P1.distance(Point(ray_position)))
                        P2 = Point(string.coords[-1][:2])
                        Pn.append(P2)
                        Pn_dist.append(P2.distance(Point(ray_position)))
                    intersection_point = Pn[Pn_dist.index(min(Pn_dist))]

                elif intersection.geom_type == 'MultiPoint':
                    print('Multipoito')

                if isinstance(closest_shape, Polygon):  # shape.geom_type == 'Polygon':

                    intersection_shape_points = list(closest_shape.exterior.coords)

                    for i in range(len(intersection_shape_points) - 1):  # ponto1 in intersection_shape_points:

                        exterior_parte = LineString(
                            [intersection_shape_points[i][:2], intersection_shape_points[i + 1][:2]])
                        dist_ref = intersection_point.distance(exterior_parte)
                        if dist_ref == intersection_point.distance(closest_shape.exterior):
                            intersection_line = exterior_parte
                            break

                    if (np.dot(dir_tx_rx, ray_direction) > 0.707) and not (
                            closest_shape in shape_ja_testados) and diffractions < max_diffractions and reflections == 0:
                        # PROCURAR QUINA
                        print('procura_quina')
                        shape_ja_testados.append(closest_shape)
                        path_difrac = tuple(path)
                        path_difrac = list(path_difrac)
                        path_difrac.append((ray_position[0], ray_position[1], 0, incidence_angle, 'reflec'))

                        for i in range(len(intersection_shape_points) - 1):  # ponto1 in intersection_shape_points:
                            p_intesec = np.array(intersection_shape_points[i][:2])
                            dir_difra = np.array(tuple(p_intesec)) - np.array(tuple(ray_position))
                            tx_quina = LineString([ray_position, p_intesec - 0.000001 * dir_difra])
                            dir_difra = dir_difra + 0.00004 * dir_difra / np.linalg.norm(dir_difra)

                            pf_difra = np.array(tuple(ray_position)) + np.array(tuple(dir_difra))
                            line = LineString([ray_position, pf_difra])

                            intersect = line.intersection(closest_shape)

                            vazio0 = False
                            if intersect.geom_type == 'Point' or intersect.is_empty:  # or intersecao.geom_type == 'GeometryCollection'
                                vazio0 = True
                                for shapedif in shapeData.geometry:
                                    intersecao_dif = tx_quina.intersection(shapedif)
                                    if not intersecao_dif.is_empty:
                                        vazio0 = False

                            if vazio0:
                                quina_rx_dir = np.array(rx_position[:2]) - np.array(tuple(p_intesec))
                                tx_rx_dir = np.array(rx_position[:2]) - np.array(tx_position[:2])
                                tx_rx_dir=tx_rx_dir/np.linalg.norm(tx_rx_dir)
                                quina_rx_dir = quina_rx_dir / np.linalg.norm(quina_rx_dir)
                                cosseno = np.dot(tx_rx_dir, quina_rx_dir)
                                if cosseno>0:
                                    print('caiu_cosseno')
                                    print(p_intesec)
                                    p_intesec2 = np.array(tuple(p_intesec)) + 0.000001 * quina_rx_dir
                                    rx_quina = LineString([rx_position[:2], p_intesec - 0.000001 * quina_rx_dir])
                                    quina_rx = LineString([p_intesec2, rx_position[:2]])
                                    intersected_segments_dif = []

                                    # intersecao_dif = GeometryCollection()
                                    intersect2 = rx_quina.intersection(closest_shape)
                                    print(intersect2)
                                    vazio = False
                                    if intersect2.geom_type == 'Point' or intersect2.is_empty:# or p_intesec.distance(intersect2) < 0.000001:
                                        vazio = True
                                        for shapedif in shapeData.geometry:
                                            intersecao_dif = quina_rx.intersection(shapedif)
                                            if not intersecao_dif.is_empty:
                                                vazio = False

                                    elif intersect2.geom_type == 'LineString':
                                        if intersect2.length<0.0000005:
                                            vazio = True
                                            for shapedif in shapeData.geometry:
                                                intersecao_dif = quina_rx.intersection(shapedif)
                                                if not intersecao_dif.is_empty:
                                                    vazio = False

                                    if vazio:
                                        print('entrou_vazio2')
                                        path_difrac.append(
                                            (p_intesec2[0], p_intesec2[1], 0, 0, 'difrac'))
                                        pnum = 1
                                        for p in ray_paths:
                                            if path_difrac == p:
                                                pnum = 0
                                        if pnum:
                                            ray_paths.append(path_difrac)
                                    else:
                                        quinas.append(p_intesec2)

                if intersection_point.distance(Point(ray_position)) >= Point(rx_position[:2]).distance(
                        Point(ray_position)):
                    print('testou hit')
                    if ray_line.distance(Point(rx_position[:2])) < max_distance:
                        print('deuhit')
                        hit_rx = True

                if reflections == 0 and reflections < max_reflections:
                    path.append((ray_position[0], ray_position[1], 0, None, 'reflec'))  # Z-component is zero for 2D
                elif reflections < max_reflections:
                    path.append((ray_position[0], ray_position[1], 0, incidence_angle,
                                 'reflec'))  # Z-component is zero for 2D
                reflected_direction, incidence_angle = reflec_raio(intersection_line, ray_position, intersection_point)
                ray_position = np.array([intersection_point.x, intersection_point.y])
                direction = reflected_direction
                reflections += 1



            else:
                if ray_line.distance(Point(rx_position[:2])) < max_distance:
                    hit_rx = True
                    path.append((rx_position[0], rx_position[1], 0, incidence_angle,
                                 'reflec'))  # Angle is zero if no reflection
                ray_position += ray_direction * 0.02
                steps += 1  # Increment the step counter

        if hit_rx:
            ray_paths.append(path)

    return ray_paths, quinas, distancia


"""# Example usage
shapefile_path = 'shapefiles/osm_urca.shp'
shapeData = load_shapefile(shapefile_path)
tx1_position = (-43.1750389,-22.9531831)#(-43.1750461, -22.9531651) #(-43.170956, -22.953926, 0)  # tx_position=(-43.1726616,-22.9536628,0)#
rx1_position =(-43.1743285,-22.9531298) #(-43.1743393, -22.9531336)  #(-43.170104, -22.953798, 0)  # rx_position = (-43.1723437,-22.9537710, 0)#
num_azimuths = 300
f = 800"""
"""ray_paths, quinas, distancia = trace_rays(shapeData, tx1_position, rx1_position, num_azimuths=num_azimuths,
                                          max_reflections=0, max_diffractions=1, diffraction_distance=1 / (3600 * 15))
"""
"""print('quinas')
#print(quinas)
tx1_position = np.array(tx1_position)
rx1_position = np.array(rx1_position)
polarizacao = 'V'
"""

def impedancia(e, u, sig, freq):
    return ((2 * np.pi * freq * 1e6 * u * 1j / (sig + 2 * np.pi * freq * 1e6 * e * 1j)) ** 0.5) / 377


def beta(e, u, sig, freq):
    raiz=(1+((sig/(2*np.pi*freq*1e6*e))**2))**0.5
    return 2 * np.pi * freq * 1e6 * ((u * e) ** 0.5) * (((1 / 2) * (raiz + 1)) ** 0.5)





def calcula_enlace(tx_position, rx_position, hg1, hg2, ray_paths, er, ersolo, sigma, sigmasolo, f, num_azimuths,polarizacao='V'):
    impe = impedancia(e0 * er, mi0, sigma, f)
    bet0 = beta(e0, mi0, 0, f)
    print(bet0)
    bet = beta(e0 * er, mi0, sigmasolo, f)
    print(bet)
    betsolo = beta(e0 * ersolo, mi0, sigmasolo, f)
    f0 = 47.7134515924  # em MHz.m
    k = f / f0

    parametros = []
    d_ref = getDistanceBetweenPointsNew(tx_position[1], tx_position[0], rx_position[1], rx_position[0])

    print('ray_paths')
    print(ray_paths)
    print(d_ref)
    cont = 0
    for gg in range(len(ray_paths)):
        path = ray_paths[gg]
        tamanho = len(path)
        d = 0
        perda = 1
        fase = 0
        if tamanho > 1:
            for i in range(1, tamanho):
                dist_ponto = np.linalg.norm(np.array(path[0][0:2]) - np.array(path[i][0:1]))
                if path[i][4] == 'reflec': # and np.linalg.norm( np.array(path[i - 1][0:2]) - np.array(path[i][0:2])) > dist_ponto * 2 * np.pi / num_azimuths:
                    print('caiu em reflec')
                    ponto_anterior = path[i - 1][0:2]
                    intersection_point, incidence_angle, ref_dif = path[i][0:2], path[i][3], path[i][4]
                    d = d + getDistanceBetweenPointsNew(ponto_anterior[1], ponto_anterior[0], intersection_point[1],
                                                        intersection_point[0])
                    if i==(tamanho-1):
                        d = d + getDistanceBetweenPointsNew(intersection_point[1], intersection_point[0], rx_position[1],
                                                            rx_position[0])

                    senfi = path[i][3]
                    x = (sigma) / (2 * np.pi * f * 1e6 * e0)
                    cos2fi = (1 - senfi ** 2)
                    if polarizacao == 'V':  # E tangente 1º casor¨r¬
                        re = (senfi - (((er - 1j * x) - cos2fi) ** 0.5)) / (senfi + (((er - 1j * x) - cos2fi) ** 0.5))
                    else:  # H tangente 2º caso r||
                        re = ((er-1j*x) * senfi - (((er - 1j * x) - cos2fi) ** 0.5)) / (
                                (er-1j*x) * senfi + (((er - 1j * x) - cos2fi) ** 0.5))

                    teta2 = np.cos(np.arcsin(((1 - incidence_angle ** 2) ** 0.5) * bet0 / bet))  # lei de snel
                    if polarizacao == 'V':  # E tangente 1º casor¨r¬
                        re = (impe * incidence_angle - teta2) / (impe * incidence_angle + teta2)
                    else:  # H tangente 2º caso r||
                        re = (incidence_angle - impe * teta2) / (incidence_angle + impe * teta2)

                    perda = perda * abs(re)
                    fase = fase + cmath.phase(re)

                elif path[i][4] == 'difrac':
                    # Achar h, d1 e d2
                    intersection_point, incidence_angle, ref_dif = path[i][0:2], path[i][3], path[i][4]
                    ponto_anterior = path[i - 1][0:2]
                    if i < tamanho - 1:
                        ponto_seguinte = path[i + 1][0:2]
                    else:
                        ponto_seguinte = list(rx_position[0:2])
                    dir1 = np.array(np.array(intersection_point) - np.array(ponto_anterior))
                    dirif = np.array(np.array(ponto_seguinte) - np.array(ponto_anterior))
                    modulo1 = np.dot(dir1, dirif) / np.linalg.norm(dirif)
                    pinter = np.array(ponto_anterior) + modulo1 * np.array(
                        np.array(intersection_point) - np.array(ponto_anterior)) / np.linalg.norm(dirif)
                    d1 = getDistanceBetweenPointsNew(ponto_anterior[1], ponto_anterior[0], pinter[1],
                                                     pinter[0])
                    d1_dif = getDistanceBetweenPointsNew(ponto_anterior[1], ponto_anterior[0], intersection_point[1],
                                                         intersection_point[0])
                    d2 = getDistanceBetweenPointsNew(ponto_seguinte[1], ponto_seguinte[0], pinter[1],
                                                     pinter[0])
                    d2_dif = getDistanceBetweenPointsNew(ponto_seguinte[1], ponto_seguinte[0], intersection_point[1],
                                                         intersection_point[0])
                    h = getDistanceBetweenPointsNew(intersection_point[1], intersection_point[0], pinter[1],
                                                    pinter[0])

                    d = getDistanceBetweenPointsNew(ponto_seguinte[1], ponto_seguinte[0], ponto_anterior[1],
                                                     ponto_anterior[0])

                    v = h * ((2 * f / 299.792458) * ((1 / d1) + (1 / d2))) ** 0.5
                    print('v')
                    print(v)
                    re = Fv(1e9)-Fv(v)
                    re = re.conjugate()
                    re = ((1+1j)/2)*re
                    perda = perda*abs(re)
                    fase = cmath.phase(re)
                    print(fase)


            d=((hg1 - hg2) ** 2 + d ** 2) ** 0.5
            d_ref2 = ((hg1 - hg2) ** 2 + d_ref ** 2) ** 0.5
            fase = fase - (d - d_ref2) * k  # verificar parson 2.19
            parametros.append([perda, fase, d])



        elif cont == 0:
            cont = 1
            perda = 1
            fase = 0
            d = d_ref
            parametros.append([perda, fase, d])

    if cont == 1:
        d = ((hg1+hg2)**2+d_ref**2)**0.5
        d_ref2=((hg1 - hg2) ** 2 + d_ref ** 2) ** 0.5
        Dphi = (d - d_ref2) *  k #k 2 * np.pi * (f / 299.792458)#4 * np.pi * hg1 * hg2 * (f / 299.792458) / d
        print(k)
        print(2 * np.pi * (f / 299.792458))
        senfi = (hg1 + hg2) / ((d ** 2 + (hg1 + hg2) ** 2) ** 0.5)
        x = sigmasolo / (2 * np.pi * f * 1e6 * e0)
        cos2fi = (1 - senfi ** 2)
        if polarizacao == 'V':  # H tangente 2º caso r||
            re = ((ersolo - 1j * x) * senfi - (ersolo - 1j * x - cos2fi) ** 0.5) / ((ersolo - 1j * x) * senfi + (ersolo - 1j * x - cos2fi) ** 0.5)
        else:  # E tangente 1º casor¨r¬
            re = (senfi - ((ersolo - 1j * x) - cos2fi) ** 0.5) / (senfi + ((ersolo - 1j * x) - cos2fi) ** 0.5)

        perda = abs(re)
        fase = cmath.phase(re) - Dphi
        parametros.append([perda, fase, d])
    print('parametros')
    print(parametros)

    ptot = 0
    for i in range(len(parametros)):
        caso=parametros[i]
        # verificar se é o caso clacular a atenuação no espaço livre antes de somar cada caso com ptot
        ptot = ptot + caso[0] * (d_ref/caso[2]) * np.exp(1j * caso[1])
    if ptot != 0:
        A = 20 * np.log10(1/abs(ptot))
    else:
        A = 0
    print('atenuação')
    print(A)
    return A


#A = calcula_enlace(tx1_position, rx1_position, h1, h2, ray_paths, er, er, sigma, sigma)

# função que recebe o caminho do shapefilie
