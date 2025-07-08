import pygame
import numpy as np
import math
import sys
from pygame import gfxdraw

# Physics Constants
A_REPULSION = 0.15    # Short-range repulsion strength
B_ATTRACTION = 0.25   # Intermediate attraction strength
C_SUPPRESSION = 0.3   # Long-range suppression strength
D_RANGE = 1.0         # Force range parameter
FORCE_SCALE = 0.015   # Force magnitude scaler
DAMPING = 0.97        # Velocity damping (energy loss)
MIN_DISTANCE = 0.1    # Minimum distance to prevent singularity
SPRING_STIFFNESS = 0.1  # Visual spring stiffness
SPRING_DAMPING = 0.95   # Visual spring damping
ENERGY_INJECTION = 3.0  # Energy added per key press

# Simulation Settings
WIDTH, HEIGHT = 900, 700
PARTICLE_RADIUS = 14
BACKGROUND_COLOR = (10, 5, 25)
PARTICLE_COLORS = [
    (255, 80, 80),    # Red quark
    (80, 255, 80),    # Green quark
    (80, 120, 255)    # Blue quark
]
SPRING_COLOR = (180, 220, 255, 180)  # Semi-transparent blue
TRAIL_LENGTH = 30
FPS = 60
COIL_COUNT = 8  # Number of coils per spring

class Spring:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        self.rest_length = np.linalg.norm(p1.pos - p2.pos)
        self.velocity = 0.0
        
    def update(self):
        # Calculate spring dynamics
        current_length = np.linalg.norm(self.p1.pos - self.p2.pos)
        stretch = current_length - self.rest_length
        
        # Spring force (Hooke's Law)
        force = SPRING_STIFFNESS * stretch
        self.velocity = (self.velocity + force) * SPRING_DAMPING
        
        # Apply forces to particles
        direction = (self.p2.pos - self.p1.pos) / (current_length + 1e-5)
        self.p1.vel += direction * self.velocity * 0.5
        self.p2.vel -= direction * self.velocity * 0.5

class Particle:
    def __init__(self, x, y, color_idx):
        self.pos = np.array([x, y], dtype=float)
        self.vel = np.array([0.0, 0.0], dtype=float)
        self.acc = np.array([0.0, 0.0], dtype=float)
        self.color_idx = color_idx
        self.trail = []
        self.connections = []
        
    def update(self, dt):
        self.vel += self.acc * dt
        self.vel *= DAMPING
        self.pos += self.vel * dt
        self.acc = np.array([0.0, 0.0])
        
        # Keep trail of positions
        self.trail.append(self.pos.copy())
        if len(self.trail) > TRAIL_LENGTH:
            self.trail.pop(0)
            
        # Boundary conditions with momentum conservation
        if self.pos[0] < PARTICLE_RADIUS:
            self.pos[0] = PARTICLE_RADIUS
            self.vel[0] = -self.vel[0] * 0.8
        elif self.pos[0] > WIDTH - PARTICLE_RADIUS:
            self.pos[0] = WIDTH - PARTICLE_RADIUS
            self.vel[0] = -self.vel[0] * 0.8
            
        if self.pos[1] < PARTICLE_RADIUS:
            self.pos[1] = PARTICLE_RADIUS
            self.vel[1] = -self.vel[1] * 0.8
        elif self.pos[1] > HEIGHT - PARTICLE_RADIUS:
            self.pos[1] = HEIGHT - PARTICLE_RADIUS
            self.vel[1] = -self.vel[1] * 0.8
            
    def add_energy(self, direction):
        """Add energy to the quark in a specific direction"""
        self.vel += direction * ENERGY_INJECTION

def strong_force(r):
    """Modified Yukawa potential force function"""
    if r < MIN_DISTANCE:
        return 30.0  # Strong repulsion at very close range
    return A_REPULSION/(r**2) - B_ATTRACTION*math.exp(-r/D_RANGE)/r - C_SUPPRESSION*math.exp(-r/0.5)

def calculate_forces(particles):
    """Compute forces between all particle pairs"""
    n = len(particles)
    for i in range(n):
        for j in range(i + 1, n):
            r_vec = particles[j].pos - particles[i].pos
            distance = np.linalg.norm(r_vec)
            
            if distance > 0:
                direction = r_vec / distance
                force_mag = strong_force(distance)
                
                # Apply equal and opposite forces
                force = force_mag * FORCE_SCALE * direction
                particles[i].acc += force
                particles[j].acc -= force

def draw_spring(surface, start, end, color, width=2, coils=COIL_COUNT):
    """Draw a spring between two points with dynamic coils"""
    # Calculate vector between points
    vector = end - start
    distance = np.linalg.norm(vector)
    
    # Skip drawing if too close
    if distance < 5:
        return
    
    # Normalize direction vector
    if distance > 0:
        direction = vector / distance
    else:
        direction = np.array([1.0, 0.0])
        
    # Calculate perpendicular vector for coil displacement
    perp = np.array([-direction[1], direction[0]])
    
    # Calculate coil parameters based on distance
    coil_count = max(3, min(20, int(coils * distance / 100)))
    segment_length = distance / coil_count
    coil_amplitude = min(8.0, segment_length * 0.3)
    
    # Create points along the spring path
    points = [start]
    
    # Add coils
    for i in range(1, coil_count):
        t = i / coil_count
        coil_phase = math.pi * 2 * i
        
        # Calculate position along the line
        base_pos = start + direction * (t * distance)
        
        # Calculate coil offset
        coil_offset = perp * math.sin(coil_phase) * coil_amplitude
        
        # Add coil point
        points.append(base_pos + coil_offset)
    
    points.append(end)
    
    # Draw the spring with anti-aliasing
    for i in range(len(points) - 1):
        pygame.draw.aaline(surface, color, 
                          points[i].astype(int), 
                          points[i+1].astype(int), width)

def draw_particle_trail(surface, particle):
    """Draw particle movement trail with smooth fading"""
    for i in range(1, len(particle.trail)):
        alpha = int(200 * i / len(particle.trail))
        radius = int(PARTICLE_RADIUS * 0.6 * i / len(particle.trail))
        color = (*PARTICLE_COLORS[particle.color_idx], alpha)
        pos = particle.trail[i].astype(int)
        
        # Draw anti-aliased circle for trail
        gfxdraw.filled_circle(surface, pos[0], pos[1], radius, color)
        gfxdraw.aacircle(surface, pos[0], pos[1], radius, color)

def draw_glow(surface, pos, color, radius):
    """Draw a glowing effect around particles"""
    for r in range(1, 4):
        alpha = 80 - r * 20
        glow_color = (*color, alpha)
        gfxdraw.filled_circle(surface, int(pos[0]), int(pos[1]), radius + r, glow_color)

def main():
    pygame.init()
    screen = pygame.display.set_mode((WIDTH, HEIGHT), pygame.SRCALPHA)
    pygame.display.set_caption("Proton Simulation - Strong Nuclear Force with Dynamic Springs")
    clock = pygame.time.Clock()
    
    # Create three particles in triangular formation
    particles = [
        Particle(WIDTH//2, HEIGHT//2 - 50, 0),  # Red quark
        Particle(WIDTH//2 - 45, HEIGHT//2 + 40, 1),  # Green quark
        Particle(WIDTH//2 + 45, HEIGHT//2 + 40, 2)   # Blue quark
    ]
    
    # Create springs between all particles
    springs = [
        Spring(particles[0], particles[1]),
        Spring(particles[1], particles[2]),
        Spring(particles[2], particles[0])
    ]
    
    # Give initial momentum
    particles[0].vel = np.array([1.5, 1.0])
    particles[1].vel = np.array([-1.0, -1.5])
    particles[2].vel = np.array([0.7, -0.8])
    
    # Font for UI
    font = pygame.font.SysFont('Arial', 16)
    title_font = pygame.font.SysFont('Arial', 24, bold=True)
    
    # Create surface for glow effects
    glow_surface = pygame.Surface((WIDTH, HEIGHT), pygame.SRCALPHA)
    
    # Main loop
    running = True
    while running:
        dt = 0.1  # Fixed time step for stability
        
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_r:  # Reset simulation
                    particles = [
                        Particle(WIDTH//2, HEIGHT//2 - 50, 0),
                        Particle(WIDTH//2 - 45, HEIGHT//2 + 40, 1),
                        Particle(WIDTH//2 + 45, HEIGHT//2 + 40, 2)
                    ]
                    springs = [
                        Spring(particles[0], particles[1]),
                        Spring(particles[1], particles[2]),
                        Spring(particles[2], particles[0])
                    ]
                elif event.key == pygame.K_UP:
                    # Add upward energy to red quark
                    particles[0].add_energy(np.array([0, -1]))
                elif event.key == pygame.K_DOWN:
                    # Add downward energy to red quark
                    particles[0].add_energy(np.array([0, 1]))
                elif event.key == pygame.K_LEFT:
                    # Add left energy to green quark
                    particles[1].add_energy(np.array([-1, 0]))
                elif event.key == pygame.K_RIGHT:
                    # Add right energy to green quark
                    particles[1].add_energy(np.array([1, 0]))
                elif event.key == pygame.K_w:
                    # Add upward energy to blue quark
                    particles[2].add_energy(np.array([0, -1]))
                elif event.key == pygame.K_s:
                    # Add downward energy to blue quark
                    particles[2].add_energy(np.array([0, 1]))
                elif event.key == pygame.K_ESCAPE:
                    running = False
        
        # Physics update
        calculate_forces(particles)
        
        # Update spring dynamics
        for spring in springs:
            spring.update()
        
        for particle in particles:
            particle.update(dt)
        
        # Rendering
        screen.fill(BACKGROUND_COLOR)
        glow_surface.fill((0, 0, 0, 0))
        
        # Draw springs
        for spring in springs:
            dist = np.linalg.norm(spring.p1.pos - spring.p2.pos)
            alpha = max(50, min(220, int(255 - dist * 1.5)))
            width = max(1, min(4, int(5 - dist/30)))
            spring_color = (*SPRING_COLOR[:3], alpha)
            draw_spring(screen, spring.p1.pos, spring.p2.pos, spring_color, width)
        
        # Draw particle trails
        for particle in particles:
            draw_particle_trail(screen, particle)
        
        # Draw particle glow effects
        for particle in particles:
            glow_color = PARTICLE_COLORS[particle.color_idx]
            draw_glow(glow_surface, particle.pos, glow_color, PARTICLE_RADIUS)
        
        # Apply glow
        screen.blit(glow_surface, (0, 0))
        
        # Draw particles
        for particle in particles:
            # Draw anti-aliased particle
            gfxdraw.filled_circle(screen, 
                                int(particle.pos[0]), 
                                int(particle.pos[1]), 
                                PARTICLE_RADIUS, 
                                PARTICLE_COLORS[particle.color_idx])
            gfxdraw.aacircle(screen, 
                           int(particle.pos[0]), 
                           int(particle.pos[1]), 
                           PARTICLE_RADIUS, 
                           PARTICLE_COLORS[particle.color_idx])
            
            # Draw inner highlight
            highlight_pos = particle.pos - np.array([PARTICLE_RADIUS*0.3, PARTICLE_RADIUS*0.3])
            gfxdraw.filled_circle(screen, 
                                int(highlight_pos[0]), 
                                int(highlight_pos[1]), 
                                int(PARTICLE_RADIUS*0.4), 
                                (255, 255, 255, 100))
        
        # Draw UI
        ui_text = [
            "STRONG NUCLEAR FORCE SIMULATION",
            "Particles: 3 'quarks' confined in a proton by dynamic springs",
            "Controls: [R] Reset  [ESC] Quit",
            f"Force Model: F(r) = {A_REPULSION}/r² - {B_ATTRACTION}·exp(-r/{D_RANGE})/r - {C_SUPPRESSION}·exp(-r/0.5)",
            "Add Energy:",
            "  Red Quark: [UP]/[DOWN] arrows",
            "  Green Quark: [LEFT]/[RIGHT] arrows",
            "  Blue Quark: [W]/[S] keys"
        ]
        
        # Draw title
        title_surface = title_font.render("QUARK CONFINEMENT IN A PROTON", True, (200, 220, 255))
        screen.blit(title_surface, (WIDTH//2 - title_surface.get_width()//2, 15))
        
        # Draw info
        for i, text in enumerate(ui_text):
            text_surface = font.render(text, True, (180, 200, 230))
            screen.blit(text_surface, (20, 60 + i*25))
        
        # Draw FPS
        fps_text = font.render(f"FPS: {int(clock.get_fps())}", True, (150, 220, 150))
        screen.blit(fps_text, (WIDTH - 120, 20))
        
        pygame.display.flip()
        clock.tick(FPS)

if __name__ == "__main__":
    main()
    pygame.quit()
    sys.exit()
